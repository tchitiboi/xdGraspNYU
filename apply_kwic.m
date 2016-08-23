function [recon_kwic] = apply_kwic(kdata,Traj,DensityComp,b1,firstRingSize,Nproj,bFilter,bSW)
%function [recon_kwic,kdata_Under,Traj_Under,DensityComp_Under,kwicmask,kwicdcf] = apply_kwic(kdata,Traj,DensityComp,b1,firstRingSize,Nproj,bFilter,bSW)

% Required parameters
%Nproj           =  89;
%firstRingSize   = 15*2; % Number of projections for first ring size
nyqFactor       = 1;%.5; % Determines radii of KWIC rings (1 ensures full Nyquist sampling for each ring)

centerP         = floor(Nproj/2)+1; % Allows asymmetric KWIC-filters
level           = 0.5; % Level for Tukey-filter (0 = no smoothing)

%% Calculate KWIC mask based on Thomas Benkert's routine
[kwicmask, kwicdcf]     = KWICnyu(size(kdata,1), Nproj, nyqFactor, firstRingSize, centerP, level);

if(bSW)
    kwicmask = ones(size(kwicmask));
    kwicdcf = ones(size(kwicdcf));
end

%kdata_orig = kdata;
%% 
clear kdata_Under Traj_Under DensityComp_Under
bPre = floor((Nproj-firstRingSize)/2);
kdata = padarray(kdata,[0 bPre 0]);
Traj = padarray(Traj,[0 bPre]);
DensityComp = padarray(DensityComp,[0 bPre]);
c=1;
bx2 = floor(Nproj/2);
% for ii=bx2+1:firstRingSize:size(kdata,2)-bx2
%     kdata_Under(:,:,:,c)=kdata(:,ii-bx2:ii+bx2,:);
%     Traj_Under(:,:,c)=Traj(:,ii-bx2:ii+bx2);
%     DensityComp_Under(:,:,c)=DensityComp(:,ii-bx2:ii+bx2);
%     c = c+1;
% end
for ii=bPre+1:firstRingSize:size(kdata,2)-Nproj+bPre
    kdata_Under(:,:,:,c)=kdata(:,ii-bPre:ii+Nproj-1-bPre,:);
    Traj_Under(:,:,c)=Traj(:,ii-bPre:ii+Nproj-1-bPre);
    DensityComp_Under(:,:,c)=DensityComp(:,ii-bPre:ii+Nproj-1-bPre);
    c = c+1;
end
clear kdata Traj DensityComp

[~,~,nc,nt] = size(kdata_Under);

%%
%kwicdcf = kwicdcf./max(kwicdcf(:));

%% kdata_kwic = kdata .* (kwicmask .* kwicdcf)
kdata_kwic = bsxfun(@times, kdata_Under, kwicmask);
[nx,ntviews,nc,nt]=size(kdata_Under);

clear kdata_Under kwicmask 

%param.E=MCNUFFT_kwic(Traj_Under,bsxfun(@times,DensityComp_Under,kwicdcf),b1,kwicmask);
E=MCNUFFT(Traj_Under,bsxfun(@times,DensityComp_Under,kwicdcf),b1);
clear Traj_Under DensityComp_Under kwicdcf b1

if(bFilter)
    kdata_kwic=double(squeeze(kdata_kwic).*repmat(kaiser(nx,20),[1,ntviews,nc,nt])); 
else
    kdata_kwic=double(squeeze(kdata_kwic)); 
end
recon_kwic = E'*kdata_kwic; 

% %%
% Weight1=0.01;
% param.TVWeight=max(abs(recon_GRASP_kwic(:)))*Weight1;% Cardiac dimension 
% param.TV = TV_Temp;% TV along Cardiac dimension 
% param.nite = 6;param.display = 1;
% %param.TVWeight=0;
% 
% %%%
% clc
% tic
% 
% for n=1:3
%     recon_GRASP_kwic = CSL1NlCg(recon_GRASP_kwic,param);
% end
% time=toc;
% time=time/60
% recon_GRASP_kwic=abs(single(recon_GRASP_kwic));
% 
% %recon_GRASP=recon_GRASP(113:end-112,113:end-112,:,:);
% 
% figure,imagescn(abs(recon_GRASP_kwic),[0 .003],[],[],3)


end