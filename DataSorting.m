function [kdata_Under2,Traj_Under2,DensityComp_Under2,Res_Signal_Under]=DataSorting(kdata,Traj,DensityComp,Res_Signal,nline,para);
%Data sorting into two dynamic dimensions, one cardiac and one respiratory

% [nx,ntviews,nc]=size(kdata);
% kdata_Under=zeros(nx,para.CardiacPhase*nline,nc,para.ntres);
% Traj_Under=zeros(nx,para.CardiacPhase*nline,para.ntres);
% DensityComp_Under=zeros(nx,para.CardiacPhase*nline,para.ntres);
% for ii=1:para.ntres
%     kdata_Under(:,:,:,ii)=kdata(:,(para.ES_index(ii)-para.Sys)*nline+1:(para.ES_index(ii)+para.Dia)*nline,:);
%     Traj_Under(:,:,ii)=Traj(:,(para.ES_index(ii)-para.Sys)*nline+1:(para.ES_index(ii)+para.Dia)*nline);
%     DensityComp_Under(:,:,ii)=DensityComp(:,(para.ES_index(ii)-para.Sys)*nline+1:(para.ES_index(ii)+para.Dia)*nline);
%     Res_Signal_Under(:,ii)=Res_Signal((para.ES_index(ii)-para.Sys)+1:(para.ES_index(ii)+para.Dia));
% end
% for ii=1:para.CardiacPhase-1
%     kdata_Under1(:,:,:,:,ii)=kdata_Under(:,(ii-1)*nline+1:ii*nline+nline,:,:);
%     Traj_Under1(:,:,:,ii)=Traj_Under(:,(ii-1)*nline+1:ii*nline+nline,:);
%     DensityComp_Under1(:,:,:,ii)=DensityComp_Under(:,(ii-1)*nline+1:ii*nline+nline,:);
% end
% Res_Signal_Under=Res_Signal_Under(1:end-1,:)';

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

Perc=4;
tres=floor(para.ntres/Perc);
for ii=1:Perc
    tmp1=kdata_Under1(:,:,:,(Perc-1)*tres+1:Perc*tres,:);
    tmp2=Traj_Under1(:,:,(Perc-1)*tres+1:Perc*tres,:);
    tmp3=DensityComp_Under1(:,:,(Perc-1)*tres+1:Perc*tres,:);
    
    tmp1=permute(tmp1,[1,2,4,3,5]);
    
    tmp1=reshape(tmp1,[nx,nline*tres,nc,para.CardiacPhase]);
    tmp2=reshape(tmp2,[nx,nline*tres,para.CardiacPhase]);
    tmp3=reshape(tmp3,[nx,nline*tres,para.CardiacPhase]);
    
    kdata_Under2(:,:,:,ii,:)=tmp1;
    Traj_Under2(:,:,ii,:)=tmp2;
    DensityComp_Under2(:,:,ii,:)=tmp3;
end