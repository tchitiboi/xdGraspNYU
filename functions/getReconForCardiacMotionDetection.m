
function [recon_Car,nt] = getReconForCardiacMotionDetection(kdata,Traj,DensityComp,b1,nline,para,bKwic);

%%%%%
bFilt=0;
%%%%%

nline_car=nline*2;
nt=floor(size(kdata,2)/nline_car);
[nx,ny,nc]=size(b1);
NN=floor(nx/3);
b1=b1((nx-NN)/2+1:end-(nx-NN)/2,(nx-NN)/2+1:end-(nx-NN)/2,:);
    
if(~bKwic)
    clear kdata_Under Traj_Under DensityComp_Under
    for ii=1:nt
        kdata_Under(:,:,:,ii)=kdata(:,(ii-1)*nline_car+1:ii*nline_car,:);
        Traj_Under(:,:,ii)=Traj(:,(ii-1)*nline_car+1:ii*nline_car);
        DensityComp_Under(:,:,ii)=DensityComp(:,(ii-1)*nline_car+1:ii*nline_car);
    end
    E=MCNUFFT(Traj_Under,DensityComp_Under,b1);
    [nx,ntviews,nc,nt]=size(kdata_Under);
    recon_Car=E'*double(kdata_Under.*repmat(kaiser(nx,20),[1,ntviews,nc,nt]));
    clear Traj_Under DensityComp_Under E
else
    [recon_Car,~,~,~,~,~] = apply_kwic(kdata,Traj,DensityComp,b1,nline_car,bFilt);
end
    
end