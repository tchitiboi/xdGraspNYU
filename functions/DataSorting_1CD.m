function [kdata_Under1,Traj_Under1,DensityComp_Under1,Res_Signal_Under]=DataSorting_1CD(kdata,Traj,DensityComp,Res_Signal,nline,para);
%Data sorting into one cardiac cycle only

[nx,ntviews,nc]=size(kdata);
kdata_Under=zeros(nx,para.CardiacPhase*nline,nc,para.ntres);
Traj_Under=zeros(nx,para.CardiacPhase*nline,para.ntres);
DensityComp_Under=zeros(nx,para.CardiacPhase*nline,para.ntres);
for ii=1:para.ntres
    kdata_Under(:,:,:,ii)=kdata(:,(para.ES_index(ii)-para.Sys)*nline+1:(para.ES_index(ii)+para.Dia)*nline,:);
    Traj_Under(:,:,ii)=Traj(:,(para.ES_index(ii)-para.Sys)*nline+1:(para.ES_index(ii)+para.Dia)*nline);
    DensityComp_Under(:,:,ii)=DensityComp(:,(para.ES_index(ii)-para.Sys)*nline+1:(para.ES_index(ii)+para.Dia)*nline);
    Res_Signal_Under(:,ii)=Res_Signal((para.ES_index(ii)-para.Sys)+1:(para.ES_index(ii)+para.Dia));
end
for ii=1:para.CardiacPhase
    kdata_Under1(:,:,:,:,ii)=kdata_Under(:,(ii-1)*nline+1:ii*nline,:,:);
    Traj_Under1(:,:,:,ii)=Traj_Under(:,(ii-1)*nline+1:ii*nline,:);
    DensityComp_Under1(:,:,:,ii)=DensityComp_Under(:,(ii-1)*nline+1:ii*nline,:);
end
Res_Signal_Under=Res_Signal_Under';

for ii=1:para.CardiacPhase;
    [~,index]=sort(Res_Signal_Under(:,ii),'descend');
    kdata_Under1(:,:,:,:,ii)=kdata_Under1(:,:,:,index,ii);
    Traj_Under1(:,:,:,ii)=Traj_Under1(:,:,index,ii);
    DensityComp_Under1(:,:,:,ii)=DensityComp_Under1(:,:,index,ii);
end

Perc=3;
kdata_Under1=kdata_Under1(:,:,:,1:round(para.ntres/Perc),:);
Traj_Under1=Traj_Under1(:,:,1:round(para.ntres/Perc),:);
DensityComp_Under1=DensityComp_Under1(:,:,1:round(para.ntres/Perc),:);

kdata_Under1=permute(kdata_Under1,[1,2,4,3,5]);

kdata_Under1=reshape(kdata_Under1,[nx,nline*round(para.ntres/Perc),nc,para.CardiacPhase]);
Traj_Under1=reshape(Traj_Under1,[nx,nline*round(para.ntres/Perc),para.CardiacPhase]);
DensityComp_Under1=reshape(DensityComp_Under1,[nx,nline*round(para.ntres/Perc),para.CardiacPhase]);
