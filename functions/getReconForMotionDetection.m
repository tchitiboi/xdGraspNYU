
function [recon,nt] = getReconForMotionDetection(kdata,Traj,DensityComp,b1,nline_res,NN,bKwic,Nproj,bFilt,bSW); %,V,S);

%%%%%
%bFilt=0;
%Nproj=nline_res;
%firstRingSize=15;
%%%%%

%nline_res=nline*8;
nt=floor(size(kdata,2)/nline_res);
[nx,ny,nc]=size(b1);
%NN=floor(nx/2);
b1=b1((nx-NN)/2+1:end-(nx-NN)/2,(nx-NN)/2+1:end-(nx-NN)/2,:);
    
if(~bKwic)
    clear kdata_Under Traj_Under DensityComp_Under
    for ii=1:nt
        kdata_Under(:,:,:,ii)=kdata(:,(ii-1)*nline_res+1:ii*nline_res,:);
        Traj_Under(:,:,ii)=Traj(:,(ii-1)*nline_res+1:ii*nline_res);
        DensityComp_Under(:,:,ii)=DensityComp(:,(ii-1)*nline_res+1:ii*nline_res);
    end
    E=MCNUFFT(Traj_Under,DensityComp_Under,b1);
    [nx,ntviews,nc,nt]=size(kdata_Under);
    recon=E'*double(kdata_Under.*repmat(kaiser(nx,20),[1,ntviews,nc,nt]));
    clear Traj_Under DensityComp_Under E
else
    [recon,~,~,~,~,~] = apply_kwic(kdata,Traj,DensityComp,b1,nline_res,Nproj,bFilt,bSW);
end
    
end