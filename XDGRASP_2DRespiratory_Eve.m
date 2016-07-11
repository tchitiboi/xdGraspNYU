clear all
clear classes
clc

%Load the files
cd ('C:\Users\tchitiboi\Desktop\testDataSeptum\Pt7')
load kdata.mat;
load Traj.mat;

[nx,ntviews,nc]=size(kdata);
coil=1:nc;%Coil elements used coil=[];coil=1:nc;
kdata=kdata(:,:,coil);
[nx,ntviews,nc]=size(kdata);

%Number of spokes in each cardiac phase
nline=15;

%Smoothing filter span
para.span=10;
para.flag=1;

%Cut is the number of spokes thrown away in the beginning of acquisition
%to achieve steady state
para.TR=TA/(ntviews+Cut)*nline;

% Total number of cardiac phases
para.nt=ntviews/nline;

%Generate time stamp for each cardiac phases
para=TimeStamp(para);

%%%Find coil element for Cardiac signal
%%% This may need some optimizations 
para.LF_H=0.6;para.HF_H=2;%%% initial heart rate range
para.LF_R=0.1;para.HF_R=0.4;%%% initial respiration rate range

%%%Calculate coil sensitivities
kdata=kdata.*repmat(sqrt(DensityComp),[1,1,nc]);
load ref.mat
ref=ref(:,:,coil);
[~,b1]=adapt_array_2d(squeeze(ref));
b1=double(b1/max(abs(b1(:))));clear ref

%Get respiratory motion signal
[Res_Signal,para]=GetRespiratoryMotionSignal_Block(kdata,Traj,DensityComp,b1,nline,para,0);
Res_Signal=Res_Signal./max(Res_Signal(:));

%Get cardiac motion signal
[Cardiac_Signal,para]=GetCardiacMotionSignal_Block(kdata,Traj,DensityComp,b1,nline,para);
Cardiac_Signal=Cardiac_Signal./max(Cardiac_Signal(:));

para=ImproveCardiacMotionSignal(Cardiac_Signal,para);

% % trying to evaluate the trajectories for the whole data set
% Nread = nx;
% Nproj=ntviews;
% rho =  linspace(0,Nread-1,Nread)' - (Nread-1)/2;
% rho = rho/Nread;
% GA = 111.246117975/180*pi;
% phi = [1:Nproj]*GA;
% traj = rho * exp(1j*phi);
% figure, plot(traj)

% % %code for 9 cardiac phases 9 resp phases
Perr=9;
Perc=9;
[kdata_Under,Traj_Under,DensityComp_Under,Res_Signal_Under]=DataSorting_Resp(kdata,Traj,DensityComp,Res_Signal,nline,para, Perc, Perr);

% [kdata_Under,Traj_Under,DensityComp_Under,Res_Signal_Under]=DataSorting_1CD(kdata,Traj,DensityComp,Res_Signal,nline,para);


%param.E=MCNUFFT_MP(Traj_Under,DensityComp_Under,b1);
% param.E=MCNUFFT(Traj_Under,DensityComp_Under,b1);
param.E=MCNUFFT_MP(Traj_Under,DensityComp_Under,b1);
%param.E=MCNUFFT(Traj_Under,DensityComp_Under,b1);

param.y=double(squeeze(kdata_Under));
% param.Res_Signal=Res_Signal_Under;
recon_GRASP=param.E'*param.y;

Weight1=0.01; %Cardiac
Weight2=0.03; %Resp
param.TVWeight=max(abs(recon_GRASP(:)))*Weight1;% Cardiac dimension 
param.L1Weight=max(abs(recon_GRASP(:)))*Weight2;% Respiration
% param.TV = TV_Temp;% TV along Cardiac dimension 
% param.TV = TV_Temp3D;% TV along Cardiac dimension 
% param.W  = TV_Temp2DRes;% TV along Respiratory dimension
% param.nite = 6;param.display = 1;
% param.TV = TV_Temp;% TV along Cardiac dimension 
param.TV = TV_Temp3D;% TV along Cardiac dimension 
param.W  = TV_Temp2DRes;% TV along Respiratory dimension
param.nite = 6;param.display = 1;


clear para Cardiac_Signal Cut DensityComp DensityComp_Under
clear Gating_Signal Gating_Signal_FFT Res_Signal Res_Signal_Under
clear TA Traj Traj_Under Weight1 Weight2 b1 kdata kdata_Under nc
clear nline ntviews nx N ans

%%%
clc
tic
for n=1:3
    recon_GRASP = CSL1NlCg(recon_GRASP,param);
end
time=toc;
time=time/60
recon_GRASP=abs(single(recon_GRASP));

[nx,ny,nt]=size(recon_GRASP);
for ii=1:nx
  for jj=1:ny
    recon_GRASP_TM(ii,jj,:)=medfilt1(recon_GRASP(ii,jj,:),5);
 end
end

figure,imagescn(abs(recon_GRASP),[0 .003],[],[],3)

figure,imagescn(abs(recon_GRASP_TM),[0 .003],[],[],3)