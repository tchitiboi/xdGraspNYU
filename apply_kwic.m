function [recon_kwic,kdata_Under,Traj_Under,DensityComp_Under,kwicmask,kwicdcf] = apply_kwic(kdata,Traj,DensityComp,b1,firstRingSize,bFilter)

% Required parameters
Nproj           = 89;
%firstRingSize   = 15*2; % Number of projections for first ring size
nyqFactor       = .5; % Determines radii of KWIC rings (1 ensures full Nyquist sampling for each ring)

centerP         = floor(Nproj/2)+1; % Allows asymmetric KWIC-filters
level           = 0.5; % Level for Tukey-filter (0 = no smoothing)

%% Calculate KWIC mask based on Thomas Benkert's routine
[kwicmask, kwicdcf]     = KWICnyu(size(kdata,1), Nproj, nyqFactor, firstRingSize, centerP, level);

%% 
clear kdata_Under Traj_Under DensityComp_Under
bx2 = floor(Nproj/2);
size(bx2+1:firstRingSize:size(kdata,2)-bx2,2);
c=1;
for ii=bx2+1:firstRingSize:size(kdata,2)-bx2
    kdata_Under(:,:,:,c)=kdata(:,ii-bx2:ii+bx2,:);
    Traj_Under(:,:,c)=Traj(:,ii-bx2:ii+bx2);
    DensityComp_Under(:,:,c)=DensityComp(:,ii-bx2:ii+bx2);
    c = c+1;
end

[~,~,nc,nt] = size(kdata_Under);

%%
kwicdcf = kwicdcf./max(kwicdcf(:));

%% kdata_kwic = kdata .* (kwicmask .* kwicdcf)
kdata_kwic = bsxfun(@times, kdata_Under, kwicmask);

param.E=MCNUFFT_kwic(Traj_Under,bsxfun(@times,DensityComp_Under,kwicdcf),b1,kwicmask);
%param.E=MCNUFFT(Traj_Under,bsxfun(@times,DensityComp_Under,kwicdcf),b1);

[nx,ntviews,nc,nt]=size(kdata_Under);
if(bFilter)
    param.y=double(squeeze(kdata_kwic).*repmat(kaiser(nx,20),[1,ntviews,nc,nt]));
else
    param.y=double(squeeze(kdata_kwic));
end
recon_kwic = param.E'*param.y; 

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