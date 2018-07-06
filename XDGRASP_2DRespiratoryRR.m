%%%%%%%%%%%%%
load_files=0;
%%%%%%%%%%%%%

addpath(genpath('functions'))
addpath(genpath('imageAnalysis'))

%cd ('/home/chitit01/NYUShared/axell01lab/labspace/Teodora/CardiacRadialData/2DRadialCardiac#SAeraTisch#F135584#M3363#D201217#T192124#Cardiac_Radial_Adult_RDD_lowEF_a/')

load kdata.mat;
load Traj.mat;
load ref.mat

%coil=[1,2,3,4,5,6,7,8,9,10];%Coil elements used coil=[];
% coil=[1:11,13:22,25:34];
% kdata=kdata(:,:,coil);

%coil compression
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
coil=1:nc;%Coil elements used coil=[];coil=1:nc;
kdata=kdata(:,:,coil);
[nx,ntviews,nc]=size(kdata);

%Number of spokes in each cardiac phase
nline=15;

% %crop data
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
b1=double(b1/max(abs(b1(:))));%clear ref

%%%%%%%%%%%%%
% ResSort=0;Teodora/meas_MID01554_FID34299_2D_Cardiac_RDD/
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
[Res_Signal, Res_Signal_Uninverted, para]=GetRespiratoryMotionSignal_BlockQuick(para,maskHeart,1,recon_Car,0);
Res_Signal=Res_Signal./max(Res_Signal(:));
Res_Signal_Uninverted = Res_Signal_Uninverted./max(Res_Signal_Uninverted(:)).*mean(Cardiac_Signal);

% para.LF_H=0.9; para.HF_H=1.7;%%% initial heart rate range
% TR=para.TR/2;
% time = TR:TR:para.nt*TR;
% F_S = 1/TR;F_X = 0:F_S/(para.nt-1):F_S;
% F_X=F_X-F_S/2;  %%% frequency after FFT of the motion signal
% if mod(para.nt,2)==0
%     F_X=F_X+F_X(para.nt/2);
% end
% % 
% FC_Index=find(F_X<para.HF_H & F_X>para.LF_H);
% FR_Index=find(F_X<para.HF_R & F_X>para.LF_R);
% 
% k=0;
% for step = -0.5:0.01:0.5
%     k=k+1;
%     new_signal(:,k) = Cardiac_Signal + Res_Signal_Uninverted'*step;
%     new_signal(:,k) = new_signal(:,k)/max(new_signal(:,k));
%     temp=abs(fftshift(fft(new_signal(:,k))));   
%     sign_FFT(:,k) = temp/max(temp(:));    
% end
% 
% Res_Peak=squeeze(sign_FFT(FR_Index,:));
% Car_Peak=squeeze(sign_FFT(FC_Index,:));
% 
% ratio_Peak = max(Car_Peak)./max(Res_Peak);
% n = find(ratio_Peak == max(ratio_Peak));
% %[m,n]=find(Signal_FFT==max(Signal_FFT(:)));
% figure, plot(new_signal(:,n))
% figure, plot(F_X,sign_FFT(:,n)),set(gca,'XLim',[-2 2]),set(gca,'YLim',[-.02 0.08])

old_cardiac_signal = Cardiac_Signal;
%Cardiac_Signal = new_signal(:,n);
para=ImproveCardiacMotionSignal(Cardiac_Signal,para);

para.nline=nline;
[cycleLabels, para] = LabelCycles(Cardiac_Signal, para, 'afib');

Perr = 6; %Perr=para.ntres;
%Perc=para.CardiacPhase;
Perc = 21;
[Res_Signal_Bins, Res_Signal_P] = getRespBins(Res_Signal_Uninverted', Perr);
labels = max(cycleLabels);
%[kdata_Under,Traj_Under,DensityComp_Under,Res_Signal_P_Under]=DataSorting_Resp_Cell(kdata,Traj,DensityComp,Res_Signal_Bins, Res_Signal_P, nline,para, Perr, Perc);
[kdata_Under,Traj_Under,DensityComp_Under,Res_Signal_P_Under]=DataSorting_Resp_Card_RR(kdata,Traj,DensityComp,Res_Signal_Bins, Res_Signal_P, cycleLabels, labels, nline,para, Perr, Perc);
%respBin = 2;
%[kdata_Under,Traj_Under,DensityComp_Under,Res_Signal_P_Under]=DataSorting_Card_RR(kdata,Traj,DensityComp,Res_Signal_Bins, Res_Signal_P, cycleLabels, labels, nline,para, Perr, respBin, Perc);
%[kdata_Under,Traj_Under,DensityComp_Under,Res_Signal_P_Under]=DataSorting_Resp_Cell_Label_varCycle(kdata,Traj,DensityComp,Res_Signal_Bins, Res_Signal_P, cycleLabels, labels, nline,para, Perr, Perc);
%[kdata_Under,Traj_Under,DensityComp_Under,Res_Signal_Under]=DataSorting_Resp(kdata,Traj,DensityComp,Res_Signal,nline,para, Perr, Perc);
% [kdata_Under,Traj_Under,DensityComp_Under,Res_Signal_Under]=DataSorting_1CD(kdata,Traj,DensityComp,Res_Signal,nline,para);

% param.E=MCNUFFT(Traj_Under,DensityComp_Under,b1);  
% kdata_Under2 = cell(9,10);
% Traj_Under2 = cell(9,10);
% DensityComp_Under2 = cell(9,10);
% Res_Signal_P_Under2 = cell(9,10);
% 
% kdata_Under2(7:9,:) = kdata_Under;
% Traj_Under2(7:9,:) = Traj_Under;
% DensityComp_Under2(7:9,:) = DensityComp_Under;

tic
param.E=MCNUFFT_MP_Cell_RR(Traj_Under,DensityComp_Under,b1);
%param.E=MCNUFFT_MP_Cell_RR_SingleResp(Traj_Under,DensityComp_Under,b1);
%param.E=MCNUFFT_MP_Cell(Traj_Under,DensityComp_Under,b1);
    
param.y=kdata_Under;
recon_GRASP=param.E'*param.y;
time=toc;
time=time/60

%parameters

% tic
% recon_Res = zeros(size(recon_GRASP));
% osf=2;
% wg=7;
% sw=16;
% 
% for i=1:Perr
%     for j=1:Perc
%         E = gpuNUFFT([real(col(Traj_Under{i,j})), imag(col(Traj_Under{i,j}))]',...
%             col(sqrt(DensityComp_Under{i,j})),osf,wg,sw,[nx,nx],b1,true);
% 
%         kdata_Under_Res = reshape(kdata_Under{i,j},[size(kdata_Under{i,j},1)*size(kdata_Under{i,j},2),nc]);
%         recon_Res(:,:,i,j) = E'*double(kdata_Under_Res);
%         
%         clear kdata_Under_Res
%         
%         recon_Res(:,:,i,j) = recon_Res(:,:,i,j).*size(DensityComp_Under{i,j},1)*pi/2/size(DensityComp_Under{i,j},2);
%     end
% end
% time=toc;
% time=time/60

%figure,imagescn(abs(recon_Res),[0 .003],[],[],3)
figure,imagescn(abs(squeeze(recon_GRASP(:,:,4,:,:))),[0 .003],[],[],3)

Weight3=0.006; 
Weight2=0.008;
Weight1=0.008; 
param.TVWeight=max(abs(recon_GRASP(:)))*Weight1;
param.L1Weight=max(abs(recon_GRASP(:)))*Weight2;
param.L1Weight1=max(abs(recon_GRASP(:)))*Weight3;
% param.TV = TV_Temp;% TV along Cardiac dimension 
%param.TV = TV_Temp3D;% TV along Cardiac dimension 
%param.W  = TV_Temp2DRes;% TV along Respiratory dimension
% param.nite = 6;param.display = 1;
% param.TV = TV_Temp;% TV along Cardiac dimension 
param.TV = TV_Temp4DCard;% TV along Cardiac dimension 
param.W  = TV_Temp4DRes;% TV along Respiratory dimension
param.W1  = TV_Temp4DRR;% TV along RR dimension

param.nite = 6;param.display = 1;

param.SGW = Res_Signal_P_Under;


% clear para Cardiac_Signal Cut DensityComp DensityComp_Under
% clear Gating_Signal Gating_Signal_FFT Res_Signal Res_Signal_Under
% clear TA Traj Traj_Under Weight1 Weight2 b1 kdata kdata_Under nc
% clear nline ntviews nx N ans

weightMat = zeros(size(Traj_Under));
maxSp = 0;
avgSp = 0;

% for card = 1:size(Traj_Under,1)
%     for label = 1:size(Traj_Under,2)
%         weightMat(card,label) = 1/size(Traj_Under{card, label},2);
%         if (maxSp < size(Traj_Under{card, label},2))
%             maxSp = size(Traj_Under{card, label},2);
%         end
%         avgSp = avgSp + size(Traj_Under{card, label},2);
%     end
% end
% avgSp = avgSp/(size(Traj_Under,1)*size(Traj_Under,2));

maxThresh = 100;

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

%%%
%clc
tic
for n=1:3
    recon_GRASP = CSL1NlCg_Cell_w_RR(recon_GRASP,param,weightMat);
    %recon_GRASP = CSL1NlCg_Cell_w(recon_GRASP,param);
end
time=toc;
time=time/60;
recon_GRASP1=abs(single(recon_GRASP));
%save recon_GRASP;
recon_GRASP_small = recon_GRASP1(101:end-100, 101:end-100,:,:,:);

figure,imagescn(abs(squeeze(recon_GRASP_small(:,:,1,:,:))),[0 .01],[],[],3)
figure,imagescn(abs(ipermute(squeeze(recon_GRASP_small(:,:,1,:,:)), [1 2 4 3])),[0 .01],[],[],3)

