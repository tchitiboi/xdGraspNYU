function [kdata_Under2,Traj_Under2,DensityComp_Under2,Res_Signal_Under]=DataSorting_Resp(kdata,Traj,DensityComp,Res_Signal,nline,para, nResp, nCard);
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

tcar=para.CardiacPhase/nCard;
if (tcar>=1)
  for ii=1:nCard
    kdata_Under1(:,:,:,:,ii)=kdata_Under(:,(ii-1)*nline*tcar+1:ii*nline*tcar,:,:);
    Traj_Under1(:,:,:,ii)=Traj_Under(:,(ii-1)*nline*tcar+1:ii*nline*tcar,:);
    DensityComp_Under1(:,:,:,ii)=DensityComp_Under(:,(ii-1)*nline*tcar+1:ii*nline*tcar,:);
  end
else
    disp('not enough spokes per image for this number of respiratory phases')
end
Res_Signal_Under=Res_Signal_Under';

for ii=1:nCard;
    [~,index]=sort(Res_Signal_Under(:,ii),'descend');
    kdata_Under1(:,:,:,:,ii)=kdata_Under1(:,:,:,index,ii);
    Traj_Under1(:,:,:,ii)=Traj_Under1(:,:,index,ii);
    DensityComp_Under1(:,:,:,ii)=DensityComp_Under1(:,:,index,ii);
end

%     for ii=1:nCard
%         tmp1=kdata_Under1(:,:,:,:,(ii-1)*tcar+1:ii*tcar);
%         tmp2=Traj_Under1(:,:,:,(ii-1)*tcar+1:ii*tcar);
%         tmp3=DensityComp_Under1(:,:,:,(ii-1)*tcar+1:ii*tcar);
%         
%         tmp1=permute(tmp1,[1,2,5,3,4]);
%         tmp2=permute(tmp2,[1,2,4,3]);
%         tmp3=permute(tmp3,[1,2,4,3]);
%         
%         tmp1=reshape(tmp1,[nx,nline*tcar,nc,para.ntres]);
%         tmp2=reshape(tmp2,[nx,nline*tcar,para.ntres]);
%         tmp3=reshape(tmp3,[nx,nline*tcar,para.ntres]);
%         
%         kdata_Under2(:,:,:,:,ii)=tmp1;
%         Traj_Under2(:,:,:,ii)=tmp2;
%         DensityComp_Under2(:,:,:,ii)=tmp3;
%     end

tres=para.ntres/nResp;
if (tres>=1)
    for ii=1:nResp
        tmp1=kdata_Under1(:,:,:,(ii-1)*tres+1:ii*tres,:);
        tmp2=Traj_Under1(:,:,(ii-1)*tres+1:ii*tres,:);
        tmp3=DensityComp_Under1(:,:,(ii-1)*tres+1:ii*tres,:);
        
        tmp1=permute(tmp1,[1,2,4,3,5]);
        
        tmp1=reshape(tmp1,[nx,size(kdata_Under1,2)*floor(tres),nc,nResp]);
        tmp2=reshape(tmp2,[nx,size(kdata_Under1,2)*floor(tres),nResp]);
        tmp3=reshape(tmp3,[nx,size(kdata_Under1,2)*floor(tres),nResp]);
        
        kdata_Under2(:,:,:,ii,:)=tmp1;
        Traj_Under2(:,:,ii,:)=tmp2;
        DensityComp_Under2(:,:,ii,:)=tmp3;
    end
else
    disp('not enough spokes per image for this number of cardiac phases')
end


end%of function