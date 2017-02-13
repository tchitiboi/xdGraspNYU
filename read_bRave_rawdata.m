%%%%%% READ IN DATA ACQUIRED WITH bRAVE %%%%%%% %%%
% author:   Rebecca Ramb
% date:     2016-07-26
% mail:     rebecca.ramb@nyumc.org
% based on Thomas Benkert's RAVE sequence, readout routine and recon
% based on Li Feng's XD-GRASP routine 
%%%%% %%%%% %%%%% %%%%% %%%%% %%%%% %%%%% %%%%% %%%


function [kdata,Traj,DensityComp,ref,TA,Cut,param] = read_bRave_rawdata(datapath,saveoutput,angle)

%%%%%%%%
Cut=480;
%%%%%%%%

if nargin < 3
    angle = 1;
end

if nargin < 2
    saveoutput=0;
end


twix_obj = mapVBVD_rr(datapath,'readheader');
d =  size(twix_obj,2);

header      = twix_obj{d}.hdr;
data        = twix_obj{d}.image;


TA=header.MeasYaps.lTotalScanTimeSec;


try
    param.noise   = twix_obj{d}.noise.unsorted();
    fprintf('Read in noise data \n')
    
catch
    fprintf('No noise data acquired \n')
    
end


try
    if strcmp(header.Meas.SoftwareVersions(end-2),'B') || strcmp(header.Meas.SoftwareVersions(end-3),'B')
        disp('Software version: VB')
        param.isVD = false;
    else
        disp('Software version: VD')
        param.isVD = true;
    end
catch
    disp('Software version could not be detected, assuming: VD')
    param.isVD = true;
end

try
    disp('Read in calibration data')
    param.calib           = twix_obj{d}.calib.unsorted();
    
    if param.isVD
        param.calibparamam1 = twix_obj{d}.calib.iceparamam(5,:); % 0/90/180/270deg scans
        param.calibparamam2 = twix_obj{d}.calib.iceparamam(6,:); % repetition
    else
        param.calibparamam1 = twix_obj{d}.calib.freeparamam(1,:); % 0/90/180/270deg scans
        param.calibparamam2 = twix_obj{d}.calib.freeparamam(2,:); % repetition
    end
catch
    disp('No calibration data acquired or error while reading in data')
    param.gradcalibration = 0;
end


%% Read the sequence params from data
param.Nread       = double(data.NCol);
param.Nproj       = double(data.NLin);
param.Ncoil       = double(data.NCha);
param.Nslc        = double(data.NSli);
param.Neco        = double(data.NEco);

%param.fftsize     = floor([param.Nread/param.osfactor, param.Nread/param.osfactor]);

%% Get k-space data
kdata = squeeze(twix_obj{d}.image());
kdata=permute(squeeze(kdata),[1,3,2]);

%% Calculate trajectory data
fprintf('Calculate trajectory \n')
%[Traj,DensityComp] = calcTrajectory(par);
%[Traj,DensityComp]=Trajectory_GoldenAngle(param.Nproj,param.Nread);
[Traj,DensityComp]=Trajectory_GoldenAngle_var(param.Nproj,param.Nread,angle);
param.angle = angle;


%% Cut the first 480 spokes that are in non-steady-state
% Cut=480;
kdata=kdata(:,Cut+1:end,:,:);
Traj=Traj(:,Cut+1:end);
DensityComp=DensityComp(:,Cut+1:end);
[param.Nread,param.Nproj,param.Ncoil]=size(kdata);

%Do some filtering that are needed to calculated coil sensitivities
filter=kaiser(param.Nread,10);
kdata1=kdata.*repmat(filter,[1,param.Nproj,param.Ncoil]);
kdata1=kdata1.*repmat(sqrt(DensityComp),[1,1,param.Ncoil]);

%Reconstruct the averaged results
if(size(Traj,2)<=6000)
    N=size(Traj,2);
else
    N=6000; %use to generate ref
end
% N=30; %use to check artifacts in each individual coil element
param.E = MCNUFFT(Traj(:,1:N),DensityComp(:,1:N),ones(param.Nread,param.Nread));
for ch=1:param.Ncoil
    ch
    ref(:,:,ch) = param.E'*double(kdata1(:,1:N,ch));
end
ref=single(ref/max(ref(:)));

%Results from individual coil element
figure,imagescn(abs(ref),[0 .5],[],[],4)

%Results combining all the coil elements
figure,imagescn(abs(sos(ref,3)),[0 1],[],[],4)


    
end