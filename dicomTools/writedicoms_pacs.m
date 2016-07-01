%% Exports dicoms to PACS
% Procedure:
% - Export single DICOM from PACS from respective patient
% - Get DICOM tags from this file
% - Overwrite specific tags such as seriesnumber etc.
% - Export dicoms
% - Push dicoms to PACS via command:
% for f in *.dcm; do storescu -to 60 -ta 60 -aet CBI_GRASP -aec NYURAD_SERIES 10.185.67.100 104 „$f“; done

%%
clear;
clc;

%% Choose measurement

% cd 'Y:\homes\benket01\Measdata\00_BREAST\20160313_DCE'
% filename = 'fatwater_MID1888_3D_GRASP_proj21_tv08';

% cd 'Y:\homes\benket01\Measdata\00_NECK\20160402'
% filename = 'fatwater_MIDThomasNeck#SSystem9#F322549#M1246#D010416#T075655#RAVE_dixon_TBN_3D_itergn4_demodfm';

cd 'Y:\homes\benket01\Measdata\00_NECK\20160430'
filename = 'fatwater_MIDThomasNeck#SSystem9#F345100#M920#D290416#T084015#RAVE_dixon_TBN_3D_proj600';

% cd 'Y:\homes\benket01\Measdata\00_NECK\20160511'
% filename = 'fatwater_MIDThomasNeck#SSystem9#F4906#M1119#D110516#T081753#RAVE_dixon_TBN_3D_proj600';

%% Load measurement
% load(filename, 'water', 'par', 'header')
load(filename, 'water', 'fat', 'par', 'header')

%% Scale images

foldername = 'F345100_dicoms_rave';
% Nomenclature: species.bin.slice
folder = [pwd '\' foldername];

% create folders
mkdir(folder)

% cast to double
water           = cast(water,'double');
fat             = cast(fat,'double');

water           = water./max(abs(water(:)))/100;
fat             = fat./max(abs(fat(:)))/100;

%%
info = dicominfo('1.3.12.2.1107.5.2.43.67041.2016042908214971360974241.dcm');

info.PatientName
header.Dicom.tPatientName

info = rmfield(info,'WindowCenter');
info = rmfield(info,'WindowWidth');

info.ProtocolName = 'RAVE_dixon_TBN';
info.SequenceName = 'RAVE3d3';

info.StudyInstanceUID = dicomuid;

info1 = info;
info2 = info;

info1.SeriesDescription = 'RAVE_dixon_water';
info2.SeriesDescription = 'RAVE_dixon_fat';

info1.SeriesInstanceUID = dicomuid;
info2.SeriesInstanceUID = dicomuid;

info1.SeriesNumber = 2000;
info2.SeriesNumber = 2001;

%% Save dicoms

for ip = 1:par.Npar
    % for ib = 1:par.Nframes
    % ib = 30;
    
    info1.InstanceNumber = ip;
    info2.InstanceNumber = ip;
    
    % img_show = abs(water(:,:,ib,ip));
    img_show1 = abs(water(:,:,ip));
    img_show2 = abs(fat(:,:,ip));
       
    img_show1 = imresize(img_show1,[size(img_show1,1)*4 size(img_show1,2)*4],'bilinear');
    img_show2 = imresize(img_show2,[size(img_show2,1)*4 size(img_show2,2)*4],'bilinear');
    
    % dicomwrite(img_show, [folder '\water.' num2str(ip, '%03d') '.' num2str(ib, '%03d') '.dcm'],info);
    dicomwrite(img_show1, [folder '\water.' num2str(ip, '%03d') '.dcm'],info1);
    dicomwrite(img_show2, [folder '\fat.' num2str(ip, '%03d') '.dcm'],info2);
    
    % end
end

fprintf('Wrote all dicoms \n')
