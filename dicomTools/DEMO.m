clear
clc

%Prepare for DICOM template
%You need to change the path to be the one you store the templates.
PathName1='C:/Users/chitit01/NYUdev/xdgrasp/Demo_DICOMGeneration/DICOM_Template_20Phase/'

load FileName.mat

%Load the matlab images
load recon_GRASP_Res16_TV02_TV01.mat
data= recon_GRASP/max(recon_GRASP(:));
[nx,ny,ntres,ntcar]=size(recon_GRASP);

%Go to the fold where you want to save your images
cd('C:/Users/chitit01/NYUdev/xdgrasp/Demo_DICOMGeneration/DICOM_Cardiac')

%save the images
Name='Cardiac';
for jj=1:ntres
    jj
    for ii=1:ntcar
        DI=FileName{(jj-1)*77+ii}; % 77 is hard coded, you do not have to change.
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
        IM=abs(data(:,:,jj,ii));
        [IM1, map] = gray2ind(IM,256);
        dicomwrite(IM1, [], [DI], metadata, 'CreateMode', 'copy');
    end
end