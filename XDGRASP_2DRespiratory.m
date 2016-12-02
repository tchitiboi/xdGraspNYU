%%%%%%%%%%%%%
load_files=0;
%%%%%%%%%%%%%

addpath(genpath('functions'))
addpath(genpath('imageAnalysis'))

% if(load_files)
%     clear all
%     clear classes
%     clc
%     %load('/Users/rambr01/Documents/MATLAB/xdgrasp/Pt67/Traj.mat')
%     %load('/Users/rambr01/Documents/MATLAB/xdgrasp/Pt67/ref.mat')
%     %load('/Users/rambr01/Documents/MATLAB/xdgrasp/Pt67/kdata.mat')
%     %load('/Users/rambr01/Documents/MATLAB/xdgrasp/Pt57/kdata.mat')
%     %load('/Users/rambr01/Documents/MATLAB/xdgrasp/Pt57/ref.mat')
%     %load('/Users/rambr01/Documents/MATLAB/xdgrasp/Pt57/Traj.mat')
%     %load('/Users/rambr01/Documents/MATLAB/xdgrasp/Pt85/kdata.mat')
%     %load('/Users/rambr01/Documents/MATLAB/xdgrasp/Pt85/ref.mat')
%     %load('/Users/rambr01/Documents/MATLAB/xdgrasp/Pt85/Traj.mat')
%     %load('/Users/rambr01/Documents/MATLAB/xdgrasp/Pt49/kdata.mat')
%     %load('/Users/rambr01/Documents/MATLAB/xdgrasp/Pt49/ref.mat')
%     %load('/Users/rambr01/Documents/MATLAB/xdgrasp/Pt49/Traj.mat')
%     load('/Users/rambr01/Documents/MATLAB/xdgrasp/Pt79/kdata.mat')
%     load('/Users/rambr01/Documents/MATLAB/xdgrasp/Pt79/ref.mat')
%     load('/Users/rambr01/Documents/MATLAB/xdgrasp/Pt79/Traj.mat')
% end
%Load the files
cd ('D:\Users\Public\cardiodata\testDataSeptum\LVdilated\Pt26')
load kdata.mat;
load Traj.mat;
load ref.mat

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
para.LF_H=0.8;para.HF_H=2;%%% initial heart rate range
para.LF_R=0.1;para.HF_R=0.4;%%% initial respiration rate range

%%%Calculate coil sensitivities
kdata=kdata.*repmat(sqrt(DensityComp),[1,1,nc]);
ref=ref(:,:,coil);
[~,b1]=adapt_array_2d(squeeze(ref));
b1=double(b1/max(abs(b1(:))));%clear ref

%%%%%%%%%%%%%
% ResSort=0;
% fov_size=floor(nx/2);
% bKwicResp=0;
% %%%%%%%%%%%%%
% % if(~bKwicResp)
%     nline_res = nline*4;NProj=0;
%     [recon_Res,nt] = getReconForMotionDetection(kdata,Traj,DensityComp,b1,nline_res,fov_size,0,NProj,0);
%     [Res_Signal,para]=GetRespiratoryMotionSignal_Block(kdata,Traj,DensityComp,b1,nline,nline_res,para,ResSort,nt,recon_Res);
%     %[Res_Signal_new,para,U,S,V]=GetRespiratoryMotionSignal_Block_SVD(kdata,Traj,DensityComp,b1,nline,para,ResSort);
%     %%%%%%%%%%%%%
% else
%     %bKwic=1;nline_res = nline*4; NProj=nline*10-1;
%     nline_res = nline*4; NProj=nline*9-1;
%     [recon_Res_kwic,nt_kwic] = getReconForMotionDetection(kdata,Traj,DensityComp,b1,nline_res,fov_size,1,NProj,0);
%     [Res_Signal_kwic,para]=GetRespiratoryMotionSignal_Block(kdata,Traj,DensityComp,b1,nline,nline_res,para,ResSort,nt_kwic,recon_Res_kwic);
%     Res_Signal = Res_Signal_kwic;
% end
% Res_Signal=Res_Signal./max(Res_Signal(:));
%%%%%%%%%%%%%


%Get cardiac motion signal
%%%%%%%%%%%%%
bWithMask=1; 
bKwicCard=0;
bSW=0;
%%%%%%%%%%%%%
fov_size=floor(nx/2.5); 
%%%%%%%%%%%%%
if(~bKwicCard && ~bSW)
    bFilt=1;
    nline_res=nline*2;NProj=0;
    [recon_Car,nt] = getReconForMotionDetection(kdata,Traj,DensityComp,b1,nline_res,fov_size,0,0,1,0);
else
    if(bKwicCard)
        nline_res = nline*2; NProj=nline*10-5;
        [recon_Car,nt_kwic] = getReconForMotionDetection(kdata,Traj,DensityComp,b1,nline_res,fov_size,1,NProj,0,0); 
    else
        bFilt=1;
        nline_res = nline*2; NProj=nline*10-5;
        [recon_Car,nt_sw] = getReconForMotionDetection(kdata,Traj,DensityComp,b1,nline_res,fov_size,1,NProj,bFilt,1); 
    end
end
if(bWithMask)
    [Cardiac_Signal,para,maskHeart]=GetCardiacMotionSignal_HeartBlock(kdata,Traj,DensityComp,b1,nline,para,recon_Car);
    Cardiac_Signal=Cardiac_Signal./max(Cardiac_Signal(:));
else
    [Cardiac_Signal,para]=GetCardiacMotionSignal_Block(kdata,Traj,DensityComp,b1,nline,para,recon_Car);
    Cardiac_Signal=Cardiac_Signal./max(Cardiac_Signal(:));
end

%Get respiratory motion signal
%[Res_Signal,para]=GetRespiratoryMotionSignal_Block(kdata,Traj,DensityComp,b1,nline,para,maskHeart,1);
[Res_Signal,para]=GetRespiratoryMotionSignal_BlockQuick(para,maskHeart,1,recon_Car,1);
Res_Signal=Res_Signal./max(Res_Signal(:));

para=ImproveCardiacMotionSignal(Cardiac_Signal,para);

% % %code for 9 cardiac phases 9 resp phases
Perr = 9; %Perr=para.ntres;
%Perc=para.CardiacPhase;
Perc = 9;
[Res_Signal_Bins, Res_Signal_P] = getRespBins(Res_Signal, Perr);
[kdata_Under,Traj_Under,DensityComp_Under,Res_Signal_P_Under]=DataSorting_Resp_Cell(kdata,Traj,DensityComp,Res_Signal_Bins, Res_Signal_P, nline,para, Perr, Perc);
%[kdata_Under,Traj_Under,DensityComp_Under,Res_Signal_Under]=DataSorting_Resp(kdata,Traj,DensityComp,Res_Signal,nline,para, Perr, Perc);
% [kdata_Under,Traj_Under,DensityComp_Under,Res_Signal_Under]=DataSorting_1CD(kdata,Traj,DensityComp,Res_Signal,nline,para);

% param.E=MCNUFFT(Traj_Under,DensityComp_Under,b1);    

tic
param.E=MCNUFFT_MP_Cell(Traj_Under,DensityComp_Under,b1);
%param.E=MCNUFFT_MP(Traj_Under,DensityComp_Under,b1);
    
param.y=kdata_Under;
recon_GRASP=param.E'*param.y;
time=toc;
time=time/60

%parameters

tic
recon_Res = zeros(size(recon_GRASP));
osf=2;
wg=7;
sw=16;

for i=1:Perr
    for j=1:Perc
        E = gpuNUFFT([real(col(Traj_Under{i,j})), imag(col(Traj_Under{i,j}))]',...
            col(sqrt(DensityComp_Under{i,j})),osf,wg,sw,[nx,nx],b1,true);

        kdata_Under_Res = reshape(kdata_Under{i,j},[size(kdata_Under{i,j},1)*size(kdata_Under{i,j},2),nc]);
        recon_Res(:,:,i,j) = E'*double(kdata_Under_Res);
        
        clear kdata_Under_Res
        
        recon_Res(:,:,i,j) = recon_Res(:,:,i,j).*size(DensityComp_Under{i,j},1)*pi/2/size(DensityComp_Under{i,j},2);
    end
end
time=toc;
time=time/60

figure,imagescn(abs(recon_Res),[0 .003],[],[],3)
figure,imagescn(abs(recon_GRASP),[0 .003],[],[],3)

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

param.SGW = Res_Signal_P_Under;


clear para Cardiac_Signal Cut DensityComp DensityComp_Under
clear Gating_Signal Gating_Signal_FFT Res_Signal Res_Signal_Under
clear TA Traj Traj_Under Weight1 Weight2 b1 kdata kdata_Under nc
clear nline ntviews nx N ans



%%%
clc
tic
for n=1:4
    recon_GRASP = CSL1NlCg_Cell_w(recon_GRASP,param);
end
time=toc;
time=time/60
recon_GRASP=abs(single(recon_GRASP));

figure,imagescn(abs(recon_GRASP),[0 .01],[],[],3)
figure,imagescn(abs(ipermute(recon_GRASP, [1 2 4 3])),[0 .01],[],[],3)
