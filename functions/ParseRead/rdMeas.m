function [mrprot] = rdMeas(varargin)
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
% 
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
