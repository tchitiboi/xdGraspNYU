load kdata.mat;
load Traj.mat;
load ref.mat

t_kdata= kdata;
kdata = squeeze(t_kdata(:,:,:,4));

% %coil compression
% [nx,ntviews,nc]=size(kdata);
% D = reshape(kdata, nx*ntviews,nc);
% [U,S,V] = svd(D, 'econ');
% ncc = 10;
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
%cycleLabels(:) = 1;

Perr = 7; 
Perc = 23;
[Res_Signal_Bins, Res_Signal_P] = getRespBins(Res_Signal_Uninverted', Perr);
labels = 1;
[kdata_Under,Traj_Under,DensityComp_Under,Res_Signal_P_Under]=DataSorting_Resp_Card_RR(kdata,Traj,DensityComp,Res_Signal_Bins, Res_Signal_P, cycleLabels, labels, nline,para, Perr, Perc);

tic
param.E=MCNUFFT_MP_Cell(Traj_Under,DensityComp_Under,b1);    
param.y=kdata_Under;
recon_GRASP=param.E'*param.y;
time=toc;
time=time/60

figure,imagescn(abs(recon_GRASP),[0 .003],[],[],3)

%parameters
Weight2=0.01;
Weight1=0.01; 
param.TVWeight=max(abs(recon_GRASP(:)))*Weight1;
param.L1Weight=max(abs(recon_GRASP(:)))*Weight2;
param.TV = TV_Temp3D;% TV along Cardiac dimension 
param.W  = TV_Temp2DRes;% TV along Respiratory dimension

param.nite = 6;param.display = 1;
param.b1 = b1;
param.SGW = Res_Signal_P_Under;

% clear para Cardiac_Signal Cut DensityComp DensityComp_Under
% clear Gating_Signal Gating_Signal_FFT Res_Signal Res_Signal_Under
% clear TA Traj Traj_Under Weight1 Weight2 b1 kdata kdata_Under nc
% clear nline ntviews nx N ans

%%%
clc
tic
for n=1:3
    recon_GRASP = CSL1NlCg_Cell_w(recon_GRASP,param);
end
time=toc;
time=time/60;
recon_GRASP=abs(single(recon_GRASP));

recon_GRASP_small = recon_GRASP(111:end-110, 111:end-110,:,:);
save recon_GRASP;
figure,imagescn(abs(recon_GRASP_small),[0 .01],[],[],3)
figure,imagescn(abs(ipermute(recon_GRASP_small, [1 2 4 3])),[0 .01],[],[],3)
