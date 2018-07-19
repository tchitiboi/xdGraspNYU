clear
clc

%Prepare for DICOM template
%You need to change the path to be the one you store the templates.
PathName1='C:\Users\Teodora\NYU\matlabDev\Demo_DICOMGeneration\DICOM_Templace_singleSlice\';

load FileName.mat

%pathData = 'Z:\axell01lab\labspace\EvePiekarski\MATLAB\Rawdata\ReconfromPt7\Pts7to20\Pt7\Allcoils'

pathData = 'Z:\axell01lab\labspace\EvePiekarski\MATLAB\FBBHProject\RAWDATA\RawdataCines\';
pathExport = 'Z:\axell01lab\labspace\EvePiekarski\MATLAB\FBBHProject\RAWDATA\RawdataDicomExport';

%Load the matlab images
cd (pathData)
myFiles = dir(fullfile(pathData,'*.mat')); 

for k = 1:length(myFiles)
  cd (pathData)
  dataset = myFiles(k).name
  name = strsplit(dataset, '.');
  aux = load (dataset);
  recon_GRASP = aux.(name{1});
  recon_GRASP = double(recon_GRASP);
  data= recon_GRASP/max(recon_GRASP(:));
  [nx,ny,ntcar]=size(recon_GRASP)
  newSubFolder = sprintf('%s/%s', pathExport, name{1});

  if ~exist(newSubFolder, 'dir')
    mkdir(newSubFolder);
  end
  
  %Go to the fold where you want to save your images
  cd(newSubFolder);
  %remove old results
  system('del /q *');

  for ii=1:ntcar
        DI=FileName{ii}; 
        metadata = dicominfo([PathName1 DI]);

        %Here, you can add dicom header
        metadata.Width=nx;
        metadata.Height=ny;
        metadata.SeriesDescription=name{1};
        metadata.Rows=nx;
        metadata.Columns=ny;
        metadata.SliceThickness=8; %slice thickness
        metadata.RepetitionTime=3.0; %TR
        metadata.EchoTime=1.5; %TE
        metadata.PatientName.FamilyName=name{1};
        metadata.PatientName.GivenName=[];
        metadata.PixelSpacing=[2 2]; % voxel size
        IM=abs(data(:,:,ii));
        [IM1, map] = gray2ind(IM,256);
        dicomwrite(IM1, [], [DI], metadata, 'CreateMode', 'copy');
  end  
end

