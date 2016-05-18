function [kdata_Under1,Traj_Under1,DensityComp_Under1,Res_Signal_Under]=DataSorting(kdata,Traj,DensityComp,Res_Signal,nline,para);
%Data sorting into two dynamic dimensions, one cardiac and one respiratory

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
for ii=1:para.CardiacPhase-1
    kdata_Under1(:,:,:,:,ii)=kdata_Under(:,(ii-1)*nline+1:ii*nline+nline,:,:);
    Traj_Under1(:,:,:,ii)=Traj_Under(:,(ii-1)*nline+1:ii*nline+nline,:);
    DensityComp_Under1(:,:,:,ii)=DensityComp_Under(:,(ii-1)*nline+1:ii*nline+nline,:);
end
Res_Signal_Under=Res_Signal_Under(1:end-1,:)';