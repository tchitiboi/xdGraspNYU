clear all
clear classes
clc

%Load the files
cd('G:\Datasetbadresults\Very bad result\kdata,ref,traj\Pt104\')
load kdata.mat;
load Traj.mat;

[nx,ntviews,nc]=size(kdata);
%coil=[1,2,3,4,5,6,7,8,9,10];%Coil elements used coil=[];
coil=1:nc;
kdata=kdata(:,:,coil);
[nx,ntviews,nc]=size(kdata);

%Number of spokes in each cardiac phase
nline=15;

%Smoothing filter span
para.span=10;

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
clb1=double(b1/max(abs(b1(:))));clear ref

%Get respiratory motion signal
[Res_Signal,para]=GetRespiratoryMotionSignal_Block(kdata,Traj,DensityComp,b1,nline,para,0);
Res_Signal=Res_Signal./max(Res_Signal(:));

%Get cardiac motion signal
[Cardiac_Signal,para]=GetCardiacMotionSignal_HeartBlock(kdata,Traj,DensityComp,b1,nline,para);
Cardiac_Signal=Cardiac_Signal./max(Cardiac_Signal(:));

para=ImproveCardiacMotionSignal(Cardiac_Signal,para);

%code for 4 respiratory phases
% [kdata_Under,Traj_Under,DensityComp_Under,Res_Signal_Under]=DataSorting(kdata,Traj,DensityComp,Res_Signal,nline,para);
% 
% % [kdata_Under,Traj_Under,DensityComp_Under,Res_Signal_Under]=DataSorting_1CD(kdata,Traj,DensityComp,Res_Signal,nline,para);
% 
% 
% %param.E=MCNUFFT_MP(Traj_Under,DensityComp_Under,b1);
% % param.E=MCNUFFT(Traj_Under,DensityComp_Under,b1);
% param.E=MCNUFFT_MP(Traj_Under,DensityComp_Under,b1);
% %param.E=MCNUFFT(Traj_Under,DensityComp_Under,b1);
% 
% param.y=double(squeeze(kdata_Under));
% % param.Res_Signal=Res_Signal_Under;
% recon_GRASP=param.E'*param.y;
% 
% Weight1=0.03;
% Weight2=0.01;
% param.TVWeight=max(abs(recon_GRASP(:)))*Weight1;% Cardiac dimension 
% param.L1Weight=max(abs(recon_GRASP(:)))*Weight2;% Respiration
% % param.TV = TV_Temp;% TV along Cardiac dimension 
% % param.TV = TV_Temp3D;% TV along Cardiac dimension 
% % param.W  = TV_Temp2DRes;% TV along Respiratory dimension
% % param.nite = 6;param.display = 1;
% % param.TV = TV_Temp;% TV along Cardiac dimension 
% param.TV = TV_Temp3D;% TV along Cardiac dimension 
% param.W  = TV_Temp2DRes;% TV along Respiratory dimension
% param.nite = 6;param.display = 1;
% 
% 
% clear para Cardiac_Signal Cut DensityComp DensityComp_Under
% clear Gating_Signal Gating_Signal_FFT Res_Signal Res_Signal_Under
% clear TA Traj Traj_Under Weight1 Weight2 b1 kdata kdata_Under nc
% clear nline ntviews nx N ans
% 
% %%%
% clc
% tic
% for n=1:3
%     recon_GRASP = CSL1NlCg(recon_GRASP,param);
% end
% time=toc;
% time=time/60
% recon_GRASP=abs(single(recon_GRASP));
% 
% [nx,ny,nt]=size(recon_GRASP);
% for ii=1:nx
%   for jj=1:ny
%     recon_GRASP_TM(ii,jj,:)=medfilt1(recon_GRASP(ii,jj,:),5);
%  end
% end
% 
% figure,imagescn(abs(recon_GRASP),[0 .003],[],[],3)
% 








% code for one respirtory phase
% [kdata_Under,Traj_Under,DensityComp_Under,Res_Signal_Under]=DataSorting(kdata,Traj,DensityComp,Res_Signal,nline,para);

[kdata_Under,Traj_Under,DensityComp_Under,Res_Signal_Under]=DataSorting_1CD(kdata,Traj,DensityComp,Res_Signal,nline,para);

mask_recon_GRASP = createCoilMasks(kdata_Under, Traj_Under, DensityComp_Under, b1);

% apply mask
b2 = b1.*mask_recon_GRASP;
param.E=MCNUFFT(Traj_Under,DensityComp_Under,b2);
param.y=double(squeeze(kdata_Under));
recon_GRASP = param.E'*param.y;
figure,imagescn(abs(recon_GRASP),[0 .001],[],[],3)


Weight1=0.03;
% Weight2=0.01;
param.TVWeight=max(abs(recon_GRASP(:)))*Weight1;% Cardiac dimension 
% param.L1Weight=max(abs(recon_GRASP(:)))*Weight2;% Respiration
param.TV = TV_Temp;% TV along Cardiac dimension 
% param.TV = TV_Temp3D;% TV along Cardiac dimension 
% param.W  = TV_Temp2D(param);% TV along Respiratory dimension
param.nite = 6;param.display = 1;

%clear para Cut DensityComp DensityComp_Under
%clear Gating_Signal Gating_Signal_FFT Res_Signal_Under
%clear TA Traj Traj_Under Weight1 Weight2 b1 kdata kdata_Under nc
%clear nline ntviews nx N ans hNewButton

%%%
clc
tic

for n=1:3
    recon_GRASP = CSL1NlCg(recon_GRASP,param);
end
time=toc;
time=time/60
recon_GRASP=abs(single(recon_GRASP));

%recon_GRASP=recon_GRASP(113:end-112,113:end-112,:,:);

figure,imagescn(abs(recon_GRASP),[0 .003],[],[],3)

