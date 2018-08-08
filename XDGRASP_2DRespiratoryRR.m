%%%%%%%%%%%%%
load_files=0;
%%%%%%%%%%%%%

addpath(genpath('functions'))
addpath(genpath('imageAnalysis'))

%cd ('/home/chitit01/NYUShared/axell01lab/labspace/Teodora/CardiacRadialData/2DRadialCardiac#SAeraTisch#F135584#M3363#D201217#T192124#Cardiac_Radial_Adult_RDD_lowEF_a/')

load kdata.mat;
load Traj.mat;
load ref.mat

% coil compression
%
% [nx,ntviews,nc]=size(kdata);
% nc
% D = reshape(kdata, nx*ntviews,nc);
% [U,S,V] = svd(D, 'econ');
% ncc = 25;
% for i = ncc+1:nc
%     S(i,i)=0;
% end
% newD = U*S*(V(1:ncc,:))';
% %newD = U*S(:,1:n)*V(:,1:n)';
% kdata = reshape(newD,nx,ntviews,ncc);

[nx,ntviews,nc]=size(kdata);
coil=1:nc;
kdata=kdata(:,:,coil);
[nx,ntviews,nc]=size(kdata);

%Number of spokes in each cardiac phase
nline=15;

% crop data (if you want to use less of the time acquired)
% crop_size = 10000;
% ntviews = crop_size;
% kdata = kdata(:,1:crop_size,:);
% Traj = Traj(:,1:crop_size);
% DensityComp = DensityComp(:,1:crop_size);
% [nx,ntviews,nc]=size(kdata);


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
para.LF_H=0.8;para.HF_H=1.7;%%% initial heart rate range
para.LF_R=0.2;para.HF_R=0.45;%%% initial respiration rate range

%%%Calculate coil sensitivities
kdata=kdata.*repmat(sqrt(DensityComp),[1,1,nc]);
ref=ref(:,:,coil);
[~,b1]=adapt_array_2d(squeeze(ref));
b1=double(b1/max(abs(b1(:))));

bWithMask=1; 
bKwicCard=0;
bSW=0;
fov_size=floor(nx/2.5); 

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
[Res_Signal, Res_Signal_Uninverted, para]=GetRespiratoryMotionSignal_BlockQuick(para,maskHeart,1,recon_Car,0);
Res_Signal=Res_Signal./max(Res_Signal(:));
Res_Signal_Uninverted = Res_Signal_Uninverted./max(Res_Signal_Uninverted(:)).*mean(Cardiac_Signal);

old_cardiac_signal = Cardiac_Signal;
%Cardiac_Signal = new_signal(:,n);
para=ImproveCardiacMotionSignal(Cardiac_Signal,para);

para.nline=nline;
[cycleLabels, para] = LabelCycles(Cardiac_Signal, para, 'afib');

Perr = 6; 
Perc = 21;
[Res_Signal_Bins, Res_Signal_P] = getRespBins(Res_Signal_Uninverted', Perr);
labels = max(cycleLabels);
[kdata_Under,Traj_Under,DensityComp_Under,Res_Signal_P_Under]=DataSorting_Resp_Card_RR(kdata,Traj,DensityComp,Res_Signal_Bins, Res_Signal_P, cycleLabels, labels, nline,para, Perr, Perc);

tic
param.E=MCNUFFT_MP_Cell_RR(Traj_Under,DensityComp_Under,b1);
param.y=kdata_Under;
recon_GRASP=param.E'*param.y;
time=toc;
time=time/60
figure,imagescn(abs(squeeze(recon_GRASP(:,:,4,:,:))),[0 .003],[],[],3)

% parameters

Weight3=0.008; 
Weight2=0.01;
Weight1=0.01; 
param.TVWeight=max(abs(recon_GRASP(:)))*Weight1;
param.L1Weight=max(abs(recon_GRASP(:)))*Weight2;
param.L1Weight1=max(abs(recon_GRASP(:)))*Weight3;
param.TV = TV_Temp4DCard;% TV along Cardiac dimension 
param.W  = TV_Temp4DRes;% TV along Respiratory dimension
param.W1  = TV_Temp4DRR;% TV along RR dimension

param.nite = 6; param.display = 1;
param.SGW = Res_Signal_P_Under;

% clear para Cardiac_Signal Cut DensityComp DensityComp_Under
% clear Gating_Signal Gating_Signal_FFT Res_Signal Res_Signal_Under
% clear TA Traj Traj_Under Weight1 Weight2 b1 kdata kdata_Under nc
% clear nline ntviews nx N ans

weightMat = zeros(size(Traj_Under));
maxSp = 0;
avgSp = 0;

maxThresh = 100;

% Weight matrix according to respiratory bin sizes
for resp = 1:size(Traj_Under,1)
  for card = 1:size(Traj_Under,2)
    for label = 1:size(Traj_Under,3)
        weightMat(resp,card,label) = 1/min(size(Traj_Under{resp,card,label},2),maxThresh);
        if (maxSp < size(Traj_Under{resp,card,label},2))
            maxSp = size(Traj_Under{resp,card,label},2);
        end
        avgSp = avgSp + size(Traj_Under{resp,card,label},2);
    end
  end
end
avgSp = avgSp/(size(Traj_Under,1)*size(Traj_Under,2)*size(Traj_Under,3));

weightMat = weightMat*maxThresh;

clc
tic
for n=1:3
    recon_GRASP = CSL1NlCg_Cell_w_RR(recon_GRASP,param,weightMat);
end
time=toc;
time=time/60;
recon_GRASP1=abs(single(recon_GRASP));
save recon_GRASP1;
recon_GRASP_small = recon_GRASP1(101:end-100, 101:end-100,:,:,:);

figure,imagescn(abs(squeeze(recon_GRASP_small(:,:,1,:,:))),[0 .01],[],[],3)
figure,imagescn(abs(ipermute(squeeze(recon_GRASP_small(:,:,1,:,:)), [1 2 4 3])),[0 .01],[],[],3)

