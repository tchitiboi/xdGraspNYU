clear all
clc

addpath(genpath('functions'))

%Select the rawdata file
path = '/home/chitit01/NYUShared/axell01lab/labspace/RadialMultiSlice/8873283/sl2/'
[file,path]=uigetfile(strcat(path,'\*.dat'),'Select Mat file:');

[kdata,Traj,DensityComp,ref,TA,Cut,param] = read_bRave_rawdata([path file]);

[nx,ntviews,nc,nz]=size(kdata);

%Generate sampling trajectory and density compensation function
[Traj,DensityComp]=Trajectory_GoldenAngle(ntviews,nx);

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
    ref(:,:,ch) = param.E'*double(kdata(:,1:N,ch,1));
end
ref=single(ref/max(ref(:)));

%Results from individual coil element
figure,imagescn(abs(ref),[0 .5],[],[],4)

%Results combining all the coil elements
figure,imagescn(abs(sos(ref,3)),[0 1],[],[],4)

%Save the files
cd (path)
save -v7.3 kdata.mat kdata
save Traj.mat Traj DensityComp Cut TA
save ref.mat ref
