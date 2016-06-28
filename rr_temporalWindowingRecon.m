function [recon_initial,recon_ite] = rr_temporalWindowingRecon(kdata,Traj,DensityComp,b1,param)

close all
nline=15;

param.TVWeight = 0;
param.L1Weight = 0;
Weight1=0.01;
param.TV = TV_Temp;% TV along temporal dimension 
param.nite = 2;

for i=1:10
    nline_res = floor(size(kdata,2)/(2^(i-1)));
    nt=floor(size(kdata,2)/nline_res);
    for ii=1:nt
        kdata_Under(:,:,:,ii)=kdata(:,(ii-1)*nline_res+1:ii*nline_res,:);
        Traj_Under(:,:,ii)=Traj(:,(ii-1)*nline_res+1:ii*nline_res);
        DensityComp_Under(:,:,ii)=DensityComp(:,(ii-1)*nline_res+1:ii*nline_res);
    end
    if i==1
        param.E=MCNUFFT(Traj_Under,DensityComp_Under,b1);
        recon_initial=param.E'*double(kdata_Under);
        param.TVWeight=max(abs(recon_initial(:)))*Weight1;% Cardiac dimension 
    else
        param.E=MCNUFFT(Traj_Under(:,:,1:2^(i-1)),DensityComp_Under(:,:,1:2^(i-1)),b1);
        x0 = repmat(recon_ite{i-1},[1 1 1 2]);
        x0 = reshape(x0,size(x0,1),size(x0,2),[]);
        param.y = double(kdata_Under(:,:,:,1:2^(i-1)));
        recon_initial = CSL1NlCg(x0,param);
    end
    recon_ite{i}=recon_initial; 
    clear kdata_Under Traj_Under DensityComp_Under
end

end