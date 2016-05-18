function [mrprot, mdh, fid] = rdMeas_full(varargin)
% read in raw fid data from Siemens meas.out/meas.asc file
% returns MrProt, YAPS buffer, MDH entries, and raw fid data
%   E. Auerbach, CMRR, 2012
% usage: [mrprot, mdh, fid] = rdMeas(outPath, ascPath, reInterpolate=true, deOversample=true, unused=0, doFlip=true, doRoFT=false)
% dependencies: parse_mrprot.m
%               parse_xprot.m
%               catstruct.m
%               c_str.m

version = '2012.06.14';

% ////////////////////////////////////////////////////
% ////// first, some options
% ////////////////////////////////////////////////////

% defaults
outPath = 'meas.out';
ascPath = 'meas.asc';

% EPI is often sampled on the ramps and requires interpolation
%   set reInterpolate if this should be done as it is read
reInterpolate = true;

% usually, data is oversampled 2X by default
% set deOversample if rdMeas.m should de-oversample the data
deOversample = true;

% this is a currently unused parameter retained for interface compatibility
unused = 0;

% set doFlip = true to reverse the appropriate EPI lines automatically
doFlip = true;

% set doRoFT = true to do RoFT (or not do inv FT after deoversampling...)
doRoFT = false;

% Updated Readout filter from Steen
filter = 1;

% read passed parameters
if (nargin >= 1)
    outPath = char(varargin{1});
    tpos = strfind(lower(outPath),'meas.out');
    if (tpos)
        ascPath = outPath;
        ascPath(tpos+5:tpos+7) = 'asc';
    end
    if (nargin >= 2)
        ascPath = char(varargin{2});
        if (nargin >= 3)
            reInterpolate = varargin{3};
            if (nargin >= 4)
                deOversample = varargin{4};
                if (nargin >= 5)
                    unused = varargin{5};
                    if (nargin >= 6)
                        doFlip = varargin{6};
                        if (nargin >= 7)
                            doRoFT = varargin{7};
                        end
                    end
                end
            end
        end
    end
end

% end options -------------------------

% first, open the .out/.dat file and check the header

fprintf('\nrdMeas.m version %s by eja\n------------------------------------------\n', version)
fprintf('Open %s\n',outPath)

fp = fopen(outPath, 'r', 'ieee-le');

% find file size
fseek(fp, 0, 'eof');
fsize = ftell(fp);
fprintf('File size: %.2f MB\n',fsize/(1024*1024));

fseek(fp, 0, 'bof');
hdrsize = fread(fp, 1, 'uint32');

% assume VB13 or newer to start with
VB13 = true;
VD11 = false;

% if the offset is 32 bytes, this is old (pre-VB13) data, so check for the
% meas.asc file--otherwise, the meas.asc equivalent data is embedded in the
% .dat

if (hdrsize == 32)
    % this must be pre-VB13 data
    % for pre-VB13, skip 32-byte header in meas.out and read in parameter data from meas.asc file (MrProt & YAPS data)
    VB13 = false;
    fprintf('This appears to be pre-VB13 meas.out/meas.asc data\n');
    fprintf('Open %s\n',ascPath);
    fp2 = fopen(ascPath, 'r');
    fprintf('Read MrProt & YAPS\n');
    mparr = fread(fp2, inf, 'uint8=>char');
    fprintf('Close meas.asc\n')
    fclose(fp2);
    mrprot = parse_mrprot(mparr);
elseif (hdrsize == 0)
    % this is probably >=VD11 data (multi-RAID)
    fprintf('This appears to be VD11 or newer multi-RAID data\n');
    VD11 = true;
    nMeas = fread(fp, 1, 'uint32');  % number of measurements in this multi-raid file (<=64)
    if (nMeas > 1)
        error('Found Multi-RAID file with more than one measurement -- this is currently not supported!');
    end
    % read the MrParcRaidFileEntry structures (152 byte length, always 64 present regardless of whether they are used)
    MeasID = zeros(nMeas,1,'uint32');
    FileID = zeros(nMeas,1,'uint32');
    MeasOffset = zeros(nMeas,1,'uint64');
    MeasLen = zeros(nMeas,1,'uint64');
    PatName = cell(nMeas);
    ProtName = cell(nMeas);
    for x=1:nMeas
        MeasID(x) = fread(fp,1,'*uint32');
        FileID(x) = fread(fp,1,'*uint32');
        MeasOffset(x) = fread(fp,1,'*uint64');
        MeasLen(x) = fread(fp,1,'*uint64');
        PatName{x} = fread(fp,64,'char');
        ProtName{x} = fread(fp,64,'char');
        fprintf('Meas %d: MeasID = %d; FileID = %d; Protocol: %s\n', x, MeasID(x), FileID(x), char(ProtName{x}));
    end
    fseek(fp,double(MeasOffset(1)),'bof'); % 8+(152*64)=9736, +504 to align to 512byte (MRPARCRAID_SECT_ALIGN)
    hdrsize = double(fread(fp,1,'uint32')) + double(MeasOffset(1));
end

if (VB13)
    % this must be >=VB13 data
    nEvp = fread(fp, 1, 'uint32');  % number of embedded evp files (?)
    % read all of the embedded evp files
    EvpName = cell(nEvp);
    EvpDat = cell(nEvp);
    MeasYapsIdx = 0;
    MeasIdx = 0;
    PhoenixIdx = 0;
    for x=1:nEvp
        EvpName{x} = read_cstr(fp);
        if (strcmp(char(EvpName{x}), 'MeasYaps')), MeasYapsIdx = x; end
        if (strcmp(char(EvpName{x}), 'Meas')), MeasIdx = x; end
        if (strcmp(char(EvpName{x}), 'Phoenix')), PhoenixIdx = x; end
        dsize = fread(fp, 1, 'uint32');
        EvpDat{x} = fread(fp, dsize, 'uint8=>char');
    end
    
    if (MeasIdx == 0) % this is VB13 data
        % for VB13, just find the meas.asc part and parse the mrprotocol
        fprintf('This appears to be VB13 meas.dat data\n');
        if (MeasYapsIdx == 0)
            error('meas.asc data not found within meas.dat!')
        end
        %fprintf('%s',char(EvpDat{MeasYapsIdx}))
        mrprot = parse_mrprot(char(EvpDat{MeasYapsIdx}));
    else % this is >=VB15 data
        % for VB15, we have to parse the XProtocol for YAPS parameters, since they are gone from
        % the text Phoenix protocol, which now only contains the bare minimum of parameters
        % needed for Phoenix.
        % we will still parse the text Phoenix protocol, however, since it
        % contains useful arrays that parse_xprot can not handle yet
        if (~VD11), fprintf('This appears to be VB15 or newer meas.dat data\n'); end
        fprintf('** Parsing Phoenix text protocol\n');
        if (PhoenixIdx == 0)
            error('Phoenix data not found within meas.dat!')
        end
        mrprot = parse_mrprot(char(EvpDat{PhoenixIdx}));
        
        fprintf('** Parsing Meas XProtocol\n');
        mrprot.XProtocol = parse_xprot(char(EvpDat{MeasIdx}'));
        
        % copy XProtocol MEAS and YAPS parameters to top level for backward compatibility
        mrprot = catstruct(mrprot.XProtocol.MEAS, mrprot.XProtocol.YAPS, mrprot);
    end
end

% then, the fastest way to do this seems to be to go through meas.out once
% first to read the mdh entries (fixed size), then allocate a chunk of memory
% and go back through again read the image data.  if memory isn't preallocated,
% this takes much, much longer.

fprintf('Reading first block of MDH entries\n');

% set fpos to start of mdh data
fseek(fp, hdrsize, 'bof');

% read mdh data until eof; store pointers for image data
idx = 0;
idx_max = 0;
rawidx = zeros(1,256);
status = 0;
loop_inc = 2560;
mdh = [];
did_estimate = false;
nCH = -1;
while (status == 0)
    idx = idx + 1;
    start_fpos = ftell(fp);
    
    if (VD11) % this is the new ScanHeader+ChannelHeader format (VDxx)
        % ulFlagsAndDMALength (uint32)
        mdh.ulDMALength(idx,1)                 = fread(fp,1,'ubit24=>uint32');
        mdh.ulFlags(idx,1)                     = fread(fp,1,'uint8');
        mdh.lMeasUID(idx,1)                    = fread(fp,1,'*int32');
        temp                                   = fread(fp,3,'*uint32');
        mdh.ulScanCounter(idx,1)               = temp(1);
        mdh.ulTimeStamp(idx,1)                 = temp(2);
        mdh.ulPMUTimeStamp(idx,1)              = temp(3);
        temp                                   = fread(fp,2,'*uint16');
        mdh.ushSystemType(idx,1)               = temp(1);
        mdh.ulPTABPosDelay(idx,1)              = temp(2);
        temp                                   = fread(fp,3,'*int32');
        mdh.lPTABPosX(idx,1)                   = temp(1);
        mdh.lPTABPosY(idx,1)                   = temp(2);
        mdh.lPTABPosZ(idx,1)                   = temp(3);
        temp                                   = fread(fp,3,'*uint32');
        mdh.ulReserved1(idx,1)                 = temp(1);
        mdh.aulEvalInfoMask(idx,:)             = temp(2:3);
        
        % check for ACQEND flag
        isACQEND = ( bitget(mdh.aulEvalInfoMask(idx,1),0+1) == 1 );
        if (isACQEND)
            idx = idx - 1;
            status = 1; % ACQEND
            break;
        end
        
        % check for SYNCDATA flag (e.g., VD11 AdjCoilCorr) - ignore these
        % lines for now since they use the mdh space in an different way
        % that we do not know how to decode
        isSYNCDATA = ( bitget(mdh.aulEvalInfoMask(idx,1),5+1) == 1 );
        if (isSYNCDATA)
            status = -3; % -3 == just skip this line
        end
        
        temp                                   = fread(fp,20,'*uint16');
        mdh.ushSamplesInScan(idx,1)            = temp(1);
        mdh.ushUsedChannels(idx,1)             = temp(2);
        %start sLoopCounter
        mdh.ushLine(idx,1)                     = temp(3);
        mdh.ushAcquisition(idx,1)              = temp(4);
        mdh.ushSlice(idx,1)                    = temp(5);
        mdh.ushPartition(idx,1)                = temp(6);
        mdh.ushEcho(idx,1)                     = temp(7);
        mdh.ushPhase(idx,1)                    = temp(8);
        mdh.ushRepetition(idx,1)               = temp(9);
        mdh.ushSet(idx,1)                      = temp(10);
        mdh.ushSeg(idx,1)                      = temp(11);
        mdh.ushIda(idx,1)                      = temp(12);
        mdh.ushIdb(idx,1)                      = temp(13);
        mdh.ushIdc(idx,1)                      = temp(14);
        mdh.ushIdd(idx,1)                      = temp(15);
        mdh.ushIde(idx,1)                      = temp(16);
        %end sLoopCounter
        %start sCutOffData
        mdh.ushPre(idx,1)                      = temp(17);
        mdh.ushPost(idx,1)                     = temp(18);
        %end sCutOffData
        mdh.ushKSpaceCentreColumn(idx,1)       = temp(19);
        mdh.ushCoilSelect(idx,1)               = temp(20);
        mdh.fReadOutOffcentre(idx,1)           = fread(fp,1,'*float32');
        mdh.ulTimeSinceLastRF(idx,1)           = fread(fp,1,'*uint32');
        temp                                   = fread(fp,2,'*uint16');
        mdh.ushKSpaceCentreLineNo(idx,1)       = temp(1);
        mdh.ushKSpaceCentrePartitionNo(idx,1)  = temp(2);
        %start sSliceData
        %start sSlicePosVec
        temp                                   = fread(fp,7,'*float32');
        mdh.flSag(idx,1)                       = temp(1);
        mdh.flCor(idx,1)                       = temp(2);
        mdh.flTra(idx,1)                       = temp(3);
        %end sSlicePosVec
        mdh.aflQuaternion(idx,:)               = temp(4:7);
        %end sSliceData
        temp                                   = fread(fp,30,'*uint16');
        mdh.aushIceProgramPara(idx,:)          = temp(1:24); % 24 of these now
        mdh.aushReservedPara(idx,:)            = temp(25:28); % 4 of these
        mdh.ushApplicationCounter(idx,1)       = temp(29);
        mdh.ushApplicationMask(idx,1)          = temp(30);
        mdh.ulCRC(idx,1)                       = fread(fp,1,'*uint32');
        
        % store the current file pointer in the index
        rawidx(idx) = ftell(fp);
        
        % allocate and read the channelheaders
        nCH_tmp = mdh.ushUsedChannels(idx,1);
        
        if (nCH_tmp > 64) % stop if file is corrupt
            status = -1;
            idx = idx - 1;
            fprintf('*** File is corrupt (nchannels = %d); stopping the read\n', nCH_tmp);
        end
        
        if ((nCH_tmp > 0) && (status == 0)) % nCH==0 -> ACQEND usually, or corrupt file
            if (nCH < 0)
                nCH = nCH_tmp;
            else
                if (nCH ~= nCH_tmp)
                    error('\nERROR: number of used channels changed from %d to %d in mdh #%d!\n', nCH, nCH_tmp, idx);
                end
            end
            
            for y=1:nCH
                if (y == 1)
                    % preallocate the sChannelHeader first time through
                    mdh.sCH_ulTypeAndChannelLength(idx,:) = zeros(nCH,1,'uint32');
                    mdh.sCH_lMeasUID(idx,:)               = zeros(nCH,1,'int32');
                    mdh.sCH_ulScanCounter(idx,:)          = zeros(nCH,1,'uint32');
                    mdh.sCH_ulReserved1(idx,:)            = zeros(nCH,1,'uint32');
                    mdh.sCH_ulSequenceTime(idx,:)         = zeros(nCH,1,'uint32');
                    mdh.sCH_ulUnused2(idx,:)              = zeros(nCH,1,'uint32');
                    mdh.sCH_ulChannelId(idx,:)            = zeros(nCH,1,'uint16');
                    mdh.sCH_ulUnused3(idx,:)              = zeros(nCH,1,'uint16');
                    mdh.sCH_ulCRC(idx,:)                  = zeros(nCH,1,'uint32');
                end
                
                % read the sChannelHeader
                mdh.sCH_ulTypeAndChannelLength(idx,y) = fread(fp,1,'*uint32');
                mdh.sCH_lMeasUID(idx,y)               = fread(fp,1,'*int32');
                temp                                  = fread(fp,4,'*uint32');
                mdh.sCH_ulScanCounter(idx,y)          = temp(1);
                mdh.sCH_ulReserved1(idx,y)            = temp(2);
                mdh.sCH_ulSequenceTime(idx,y)         = temp(3);
                mdh.sCH_ulUnused2(idx,y)              = temp(4);
                temp                                  = fread(fp,2,'*uint16');
                mdh.sCH_ulChannelId(idx,y)            = temp(1);
                mdh.sCH_ulUnused3(idx,y)              = temp(2);
                mdh.sCH_ulCRC(idx,y)                  = fread(fp,1,'*uint32');
                
                % skip to next channelheader
                fseek(fp, double(mdh.ushSamplesInScan(idx))*2*4, 'cof');
            end
        else
            % skip to the next line
            fseek(fp, double(mdh.ushSamplesInScan(idx))*2*4, 'cof');
        end
    else
        % this is the old MDH header (VAxx, VBxx)
        mdh.ulDMALength(idx,1)                 = fread(fp,1,'*uint32');
        mdh.lMeasUID(idx,1)                    = fread(fp,1,'*int32');
        temp                                   = fread(fp,5,'*uint32');
        mdh.ulScanCounter(idx,1)               = temp(1);
        mdh.ulTimeStamp(idx,1)                 = temp(2);
        mdh.ulPMUTimeStamp(idx,1)              = temp(3);
        mdh.aulEvalInfoMask(idx,:)             = temp(4:5);

        % check for ACQEND flag
        isACQEND = ( bitget(mdh.aulEvalInfoMask(idx,1),0+1) == 1 );
        if (isACQEND)
            idx = idx - 1;
            status = 1; % ACQEND
            break;
        end

        temp                                   = fread(fp,20,'*uint16');
        mdh.ushSamplesInScan(idx,1)            = temp(1);
        mdh.ushUsedChannels(idx,1)             = temp(2);
        %start sLoopCounter
        mdh.ushLine(idx,1)                     = temp(3);
        mdh.ushAcquisition(idx,1)              = temp(4);
        mdh.ushSlice(idx,1)                    = temp(5);
        mdh.ushPartition(idx,1)                = temp(6);
        mdh.ushEcho(idx,1)                     = temp(7);
        mdh.ushPhase(idx,1)                    = temp(8);
        mdh.ushRepetition(idx,1)               = temp(9);
        mdh.ushSet(idx,1)                      = temp(10);
        mdh.ushSeg(idx,1)                      = temp(11);
        mdh.ushIda(idx,1)                      = temp(12);
        mdh.ushIdb(idx,1)                      = temp(13);
        mdh.ushIdc(idx,1)                      = temp(14);
        mdh.ushIdd(idx,1)                      = temp(15);
        mdh.ushIde(idx,1)                      = temp(16);
        %start sCutOffData
        mdh.ushPre(idx,1)                      = temp(17);
        mdh.ushPost(idx,1)                     = temp(18);
        mdh.ushKSpaceCentreColumn(idx,1)       = temp(19);
        mdh.ushDummy(idx,1)                    = temp(20);
        mdh.fReadOutOffcentre(idx,1)           = fread(fp,1,'*float32');
        mdh.ulTimeSinceLastRF(idx,1)           = fread(fp,1,'*uint32');
        temp                                   = fread(fp,10,'*uint16');
        mdh.ushKSpaceCentreLineNo(idx,1)       = temp(1);
        mdh.ushKSpaceCentrePartitionNo(idx,1)  = temp(2);
        mdh.aushIceProgramPara(idx,:)          = temp(3:6);
        mdh.aushFreePara(idx,:)                = temp(7:10);
        %start sSliceData
        %start sSlicePosVec
        temp                                   = fread(fp,7,'*float32');
        mdh.flSag(idx,1)                       = temp(1);
        mdh.flCor(idx,1)                       = temp(2);
        mdh.flTra(idx,1)                       = temp(3);
        mdh.aflQuaternion(idx,:)               = temp(4:7);
        if (VB13)
            mdh.ulChannelId(idx,1)             = fread(fp,1,'uint16=>uint32');  % now actually ushChannelId
            mdh.ushPTABPosNeg(idx,1)           = fread(fp,1,'*uint16');
        else
            mdh.ulChannelId(idx,1)             = fread(fp,1,'*uint32');
        end
        
        % store the current file pointer in the index, and skip to the next mdh
        rawidx(idx) = ftell(fp);
        fseek(fp, double(mdh.ushSamplesInScan(idx))*2*4, 'cof');
    end
    
    % first time through, estimate number of entries and preallocate memory
    if ((~did_estimate) && (idx >= 255) && (status == 0))
        idx_max = idx * ceil(fsize/(rawidx(idx)-hdrsize));
        loop_inc = floor(idx_max/100);
        if (loop_inc > 1024), loop_inc = 1024; end
        fprintf('Estimating %d & preallocating %d MDH entries\n',idx_max,idx*ceil(idx_max/idx));
        fnames = fieldnames(mdh);
        for x=1:size(fnames,1)
            fnamestr = char(fnames(x));
            mdh.(fnamestr) = repmat(mdh.(fnamestr), ceil(idx_max/idx), 1);
        end
        rawidx = repmat(rawidx, 1, ceil(idx_max/idx));
        did_estimate = true;
        fprintf('Reading MDH entries: %9d (%3d%%)',idx,round(100*(rawidx(idx)/fsize)));
    end
    
    % update waitbar periodically
    if (mod(idx,loop_inc) == 0), fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b%9d (%3d%%)',idx,round(100*(rawidx(idx)/fsize))); end
    %fprintf('\nDEBUG: idx = %9d ; rawidx(idx) = %15d ; fsize = %15d',idx,rawidx(idx),fsize);
    
    % check ending conditions
    if (ftell(fp) >= fsize)
        status = -1;
    end
    
    if (status == -1)
        fprintf('\n\n------------------------------------------------------------------\n');
        fprintf('WARNING: Found premature EOF -- this file is corrupt or truncated!\n');
        fprintf('------------------------------------------------------------------\n\n');
        idx = idx - 1;
    elseif (status == -3)
        % just discard this line, but continue through loop
        fseek(fp, double(mdh.ulDMALength(idx,1)) + start_fpos, 'bof');
        idx = idx - 1;
        status = 0;
    end
    
    %debug
    %fprintf('idx: %d, fpos: %d\n', idx, ftell(fp));
    %if (idx>50), error('stop'); end
end
%fprintf('\nDEBUG: done reading MDH\n');
if ((idx > 255) && (status > 0)), fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b%9d (%3d%%)\n',idx,100); end

% truncate mdh array if our estimate was high
if (idx < idx_max)
    fprintf('Estimate was high -- truncating MDH array\n');
    for x=1:size(fnames,1)
        fnamestr = char(fnames(x));
        mdh.(fnamestr) = mdh.(fnamestr)(1:idx,:);
    end
end

% store k-space center line in mrprot for later use
mrprot.eja_ushKSpaceCentreLineNo = max(mdh.ushKSpaceCentreLineNo);

% now that we know how many fids we have, and how big they are,
% we can allocate memory for them
if (deOversample == true)
    OSfactor = mrprot.flReadoutOSFactor;
else
    OSfactor = 1.0;
end

% for VDxx with the channelheader concept, allocate a 3D fid array; otherwise, 2D
if (~VD11)
    nCH = 1;
end

% find the maximum number of samples for actual image data
isACQEND     = (bitget(mdh.aulEvalInfoMask(:,1),1) == 1);
isRTFEEDBACK = ( (bitget(mdh.aulEvalInfoMask(:,1),2) == 1) & (bitget(mdh.aulEvalInfoMask(:,1),22) ~= 1) ); % not PHASCOR
isHPFEEDBACK = (bitget(mdh.aulEvalInfoMask(:,1),3) == 1);
if (nnz(isRTFEEDBACK)), fprintf('Found RTFEEDBACK lines -- skipping these\n'); end
if (nnz(isHPFEEDBACK)), fprintf('Found HPFEEDBACK lines -- skipping these\n'); end
maxSamplesInScan = max(mdh.ushSamplesInScan(~isACQEND & ~isRTFEEDBACK & ~isHPFEEDBACK));

fprintf('Allocating memory for raw data (%d x %d x %d complex single array)\n',maxSamplesInScan/OSfactor,idx,nCH);
fid = complex(zeros(maxSamplesInScan/OSfactor,idx,nCH,'single'));

tmp = whos('fid');
fprintf('Raw data memory successfully allocated (%.1f MB)\n',tmp.bytes/1048576);

% precompute some of the interpolation parameters if necessary
if (reInterpolate == true)
    % check if we really need to do this
    if (mrprot.alRegridMode(1) == 2)                     % 2 = REGRID_TRAPEZOIDAL
        % build the trapezoidal waveform
        rotrap = ones(1,mrprot.alRegridRampupTime(1)+mrprot.alRegridFlattopTime(1)+mrprot.alRegridRampdownTime(1),'single');
        roramp = single(0:1/(mrprot.alRegridRampupTime(1)-1):1);
        rotrap(1:mrprot.alRegridRampupTime(1))= roramp;
        rotrap(mrprot.alRegridRampupTime(1)+mrprot.alRegridFlattopTime(1)+1:end) = fliplr(roramp);
        
        % cut off the unused parts
        rotrap = rotrap(mrprot.alRegridDelaySamplesTime(1)+1:size(rotrap,2));
        rotrap = rotrap(1:floor(mrprot.aflRegridADCDuration(1))); % eja: floor added for VD11
        
        % integrate
        trapint = zeros(size(rotrap,2),1,'single');
        for z=1:size(rotrap,2)
            trapint(z) = sum(rotrap(1:z));
        end
        
        % assemble the desired k-space trajectory
        % add a point on the beginning and end since otherwise the
        %     interp1 function goes wacko
        destTraj = single(0:sum(rotrap)/(mrprot.alRegridDestSamples(1)+1):sum(rotrap));
    end
end

if ( (reInterpolate == true) && (mrprot.alRegridMode(1) >= 2) )
    str1 = ' and interpolating';
elseif (OSfactor > 1.0)
    str1 = ' and deoversampling';    
else
    str1 = '';
end
if (doRoFT)
    str2 = ' w/ forward FFT';
else
    str2 = '';
end
fprintf('Reading%s image raw data%s:   0%%', str1, str2);

% loop through the file again and read in the image data
for x=1:idx
    % set the position in the file
    %fprintf('\nDEBUG: seeking %d',rawidx(x));
    fseek(fp, rawidx(x), 'bof');
    
    for cCH=1:nCH
        % for VD11, skip the 32-byte channelheader (already read)
        if (VD11), fseek(fp,32,'cof'); end
        
        % read in the complex fid data
        nSamples = double(mdh.ushSamplesInScan(x));
        realIdx = 1:2:nSamples*2;
        imagIdx = realIdx + 1;
        [trc, fcnt] = fread(fp, nSamples*2, '*float32');
        if (fcnt ~= nSamples*2)
            fprintf('\nERROR: reading trace %d/%d failed! Expected %d samples, read %d! Will attempt to continue...\n', x, idx, nSamples*2, fcnt);
            fcnt = nSamples*2; %#ok<NASGU>
            trc = zeros(1,trc);
        end
        ctrc = complex(trc(realIdx),trc(imagIdx));
        
        % reverse EPI lines if desired
        if (doFlip == true)
            %if (findstr(dumpEvalInfoMask(mdh.aulEvalInfoMask(x,1)),'MDH_REFLECT') > 0)
            if (bitget(mdh.aulEvalInfoMask(x,1),25))
                ctrc = rot90(ctrc,2);
            end
        end
        
        % interpolate the data if necessary
        if ((reInterpolate == true) && (nSamples ~= 0))
            % check if we really need to do this
            if (mrprot.alRegridMode(1) == 2)                     % 2 = REGRID_TRAPEZOIDAL
                %if (findstr(dumpEvalInfoMask(mdh.aulEvalInfoMask(x,1)),'MDH_ACQEND') > 0)
                if (bitget(mdh.aulEvalInfoMask(x,1),1))
                    % skip postprocessing on this one
                else
                    % ***** use built-in interpolation functions
                    % assemble the actual k-space trajectory for this line
                    actualDwell = mrprot.aflRegridADCDuration(1)/nSamples;
                    zidx = 1;
                    actTraj = zeros(round(mrprot.aflRegridADCDuration(1)/actualDwell),1,'single'); % Steen added round()
                    for z = 1:actualDwell:mrprot.aflRegridADCDuration(1)
                        actTraj(zidx) = trapint(round(z));
                        zidx = zidx + 1;
                    end
                    
                    % Steen modifications
                    % interpolate  Modified on 2/07/2012 to account for
                    % the difference in sampling rate from the A/D
                    % converter
                    actTraj = interp1(1:size(trapint,1),trapint,1:actualDwell:mrprot.aflRegridADCDuration(1),'linear')';
                    if (filter == 1)
                        filterf = (cumsum(destTraj(2:end-1)) - cumsum(actTraj'));
                        filterf = filterf / sqrt(1/size(ctrc,1) * sum( abs(filterf).^2) );
                        ctrc = interp1(actTraj,ctrc.*filterf',destTraj,'linear');
                    else
                        % straight linear interpolation (old way)
                        ctrc = interp1(actTraj,ctrc,destTraj,'linear');
                    end
                    
                    ctrc = ctrc(2:end-1);
                    ctrc(isnan(ctrc)) = 0;
                    ctrc = ctrc';
                end
            end
        end
        
        if (OSfactor > 1.0)
            % de-oversample the data if desired
            fttrc = fft(ctrc);
            stpos = nSamples/(OSfactor*2.0);
            fttrc = [fttrc(1:stpos); fttrc(end-stpos+1:end)];
            if (doRoFT)
                ctrc = fftshift(fttrc);
            else
                ctrc = ifft(fttrc);
            end
        elseif (doRoFT)
            ctrc = fftshift(fft(ctrc));
        end
        
        if ( (isACQEND(x)) || (isRTFEEDBACK(x)) || (isHPFEEDBACK(x)) )
            % - don't save the data from the acqend line, since sometimes it has
            %   a large (>2000) number of samples and will inflate the entire
            %   data array unnecessarily
            % - also ignore feedback lines since they have different numbers of
            %    samples
            % - set SamplesInScan to 0 so as not to confuse the sort routine
            mdh.ushSamplesInScan(x) = 0;
        else
            fid(1:nSamples/OSfactor,x,cCH) = ctrc;
        end
    end
    
    % update waitbar periodically
    if (mod(x,loop_inc) == 0), fprintf('\b\b\b\b%3d%%',ceil(100*x/idx)); end
end
fprintf('\b\b\b\b%3d%%',100);

% clean up and finish
fprintf('\nClose meas.out\n');
fclose(fp);
fprintf('Success!\n\n');





%--------------------------------------------------------------------------
function outstr = read_cstr(fp)
% read null-terminated variable-length string from file

outstr = char(zeros(1,1000));
inchar = char(1);

idx = 1;
while (inchar ~= char(0))
    inchar = fread(fp, 1, 'uint8=>char');
    outstr(idx) = inchar;
    idx = idx + 1;
end

outstr = c_str(outstr);
