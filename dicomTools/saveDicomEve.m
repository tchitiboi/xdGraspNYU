%clear
%clc

%Prepare for DICOM template
%You need to change the path to be the one you store the templates.
PathName1='C:\Users\Teodora\NYU\matlabDev\Demo_DICOMGeneration\DICOM_Templace_singleSlice\'

load FileName.mat

pathData = 'Z:\axell01lab\labspace\Teodora\CardiacRadialData\6431788_a\arrhythmia_2\'

%pathData = 'Z:\axell01lab\labspace\EvePiekarski\MATLAB\FBBHProject\RAWDATA\RawdataDicomExport'

%Load the matlab images
cd (pathData)
%load RefCine1.mat
%recon_GRASP_crop = RefCine1;
data= recon_GRASP_crop/max(recon_GRASP_crop(:));
[nx,ny,ntcar]=size(recon_GRASP_crop)



newSubFolder = sprintf('%s/%s', pathData, 'dicom');

if ~exist(newSubFolder, 'dir')
  mkdir(newSubFolder);
end

%Go to the fold where you want to save your images
cd(newSubFolder);
%remove old results
system('del /q *');

%save the images
Name='CardiacTest';
% for jj=1:min(ntres,20)
%     jj
    for ii=1:ntcar
        DI=FileName{ii}; 
        metadata = dicominfo([PathName1 DI]);

        %Here, you can add dicom header
        metadata.Width=nx;
        metadata.Height=ny;
        metadata.SeriesDescription=Name;
        metadata.Rows=nx;
        metadata.Columns=ny;
        metadata.SliceThickness=8; %slice thickness
        metadata.RepetitionTime=3.0; %TR
        metadata.EchoTime=1.5; %TE
        metadata.PatientName.FamilyName=Name;
        metadata.PatientName.GivenName=[];
        metadata.PixelSpacing=[2 2]; % voxel size
        IM=abs(data(:,:,ii));
        [IM1, map] = gray2ind(IM,256);
        dicomwrite(IM1, [], [DI], metadata, 'CreateMode', 'copy');
    end
%end