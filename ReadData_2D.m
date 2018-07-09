clear all
clc

addpath(genpath('functions'))

%Select the rawdata file
%cd('/srv/mriscan/Archive/yarra_rds/Mobile2/RadialCardiac')
%cd('Z:\Rawdata')
% cd('/Users/Rebecca/Documents/MATLAB/xdGrasp_demo/Rawdata')
path = 'Z:\axell01lab\labspace\Teodora\CardiacRadialData\2DRadialCardiac#SAeraTisch#F215516#M5589#D300418#T212938#Cardiac_Radial_Adult_RDD_a\'
[file,path]=uigetfile(strcat(path,'\*.dat'),'Select Mat file:');
%mrprot = rdMeas_full([path file]);
[kdata,Traj,DensityComp,ref,TA,Cut,param] = read_bRave_rawdata([path file]);

%Total acquisition time (seconds)
%TA=mrprot.lTotalScanTimeSec;

%Read the rawdata
%[image_obj,MDH] = mapVBVD2014([path file]);
%kdata = image_obj.image{''};
%kdata=single(kdata);
%clear image_obj file path
%kdata=permute(squeeze(kdata),[1,3,2]);
[nx,ntviews,nc]=size(kdata);

%Generate sampling trajectory and density compensation function
[Traj,DensityComp]=Trajectory_GoldenAngle(ntviews,nx);

%Cut the first 480 spokes that are in non-steady-state
Cut=480;
kdata=kdata(:,Cut+1:end,:,:);
Traj=Traj(:,Cut+1:end);
DensityComp=DensityComp(:,Cut+1:end);
[nx,ntviews,nc]=size(kdata);

%Do some filtering that are needed to calculated coil sensitivities
filter=kaiser(nx,10);
kdata1=kdata.*repmat(filter,[1,ntviews,nc]);
kdata1=kdata1.*repmat(sqrt(DensityComp),[1,1,nc]);

%Reconstruct the averaged results
if(size(Traj,2)<=6000)
    N=size(Traj,2);
else
    N=6000; %use to generate ref
end
% N=30; %use to check artifacts in each individual coil element
param.E = MCNUFFT(Traj(:,1:N),DensityComp(:,1:N),ones(nx,nx));
for ch=1:nc
    ch
    ref(:,:,ch) = param.E'*double(kdata1(:,1:N,ch));
end
ref=single(ref/max(ref(:)));

%Results from individual coil element
figure,imagescn(abs(ref),[0 .5],[],[],4)

%Results combining all the coil elements
figure,imagescn(abs(sos(ref,3)),[0 1],[],[],4)

% Coil=[1,2,4,5,9,10,11:15,18];
% ref=ref(:,:,Coil);
% kdata=kdata(:,:,Coil);

%Save the files
cd (path)
save -v7.3 kdata.mat kdata
save Traj.mat Traj DensityComp Cut TA
save ref.mat ref
