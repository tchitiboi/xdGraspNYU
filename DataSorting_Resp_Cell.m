function [kdata_Under2,Traj_Under2,DensityComp_Under2,Res_Signal_Under]=DataSorting_Resp_Cell(kdata,Traj,DensityComp,Res_Signal,nline,para, nResp, nCard);
%Data sorting into two dynamic dimensions, one cardiac and one respiratory

[nx,ntviews,nc]=size(kdata);
kdata_Under=zeros(nx,para.CardiacPhase*nline,nc,para.ntres);
Traj_Under=zeros(nx,para.CardiacPhase*nline,para.ntres);
DensityComp_Under=zeros(nx,para.CardiacPhase*nline,para.ntres);

Res_Signal_Long = interp1( linspace(0,1,length(Res_Signal)), Res_Signal, linspace(0,1,ntviews), 'nearest');

%Separate cardiac cycles
for ii=1:para.ntres
    kdata_Under(:,:,:,ii)=kdata(:,(para.ES_index(ii)-para.Sys)*nline+1:(para.ES_index(ii)+para.Dia)*nline,:);
    Traj_Under(:,:,ii)=Traj(:,(para.ES_index(ii)-para.Sys)*nline+1:(para.ES_index(ii)+para.Dia)*nline);
    DensityComp_Under(:,:,ii)=DensityComp(:,(para.ES_index(ii)-para.Sys)*nline+1:(para.ES_index(ii)+para.Dia)*nline);
    Res_Signal_Under(:,ii)=Res_Signal_Long((para.ES_index(ii)-para.Sys)*nline+1:(para.ES_index(ii)+para.Dia)*nline);
end

%collapse cardiac phases
tcar=para.CardiacPhase/nCard;
if (tcar>=1)
  for ii=1:nCard
    kdata_Under1(:,:,:,:,ii) = kdata_Under(:,(ii-1)*nline*tcar+1:ii*nline*tcar,:,:);
    Traj_Under1(:,:,:,ii) = Traj_Under(:,(ii-1)*nline*tcar+1:ii*nline*tcar,:);
    DensityComp_Under1(:,:,:,ii) = DensityComp_Under(:,(ii-1)*nline*tcar+1:ii*nline*tcar,:);
    Res_Signal_Under1(:,:,ii) = Res_Signal_Under((ii-1)*nline*tcar+1:ii*nline*tcar,:);
  end
else
    disp('not enough spokes per image for this number of cardiac phases')
end
%Res_Signal_Under=Res_Signal_Under';
% 
% for ii=1:nCard
%     [~,index]=sort(Res_Signal_Under(:,ii),'descend');
%     kdata_Under1(:,:,:,:,ii)=kdata_Under1(:,:,:,index,ii);
%     Traj_Under1(:,:,:,ii)=Traj_Under1(:,:,index,ii);
%     DensityComp_Under1(:,:,:,ii)=DensityComp_Under1(:,:,index,ii);
% end

%initializa cell array
kdata_Under2 = cell(nResp,nCard);
Traj_Under2 = cell(nResp,nCard);
DensityComp_Under2 = cell(nResp,nCard);

tres=para.ntres/nResp;
if (tres>=1)  
    for resp = 1:nResp
        for card = 1:nCard
            r = resp;
            [ind_a, ind_b] = find(Res_Signal_Under1(:,:,card) == r);
            
            while (length(ind_a)==0 && r>0)
                r = r-1;
                [ind_a, ind_b] = find(Res_Signal_Under1(:,:,card) == r);
            end            
                
            for i = 1:length(ind_a)
                tmp1(:,i,:)=squeeze(kdata_Under1(:,ind_a(i),:,ind_b(i),card));
                tmp2(:,i)=squeeze(Traj_Under1(:,ind_a(i),ind_b(i),card));
                tmp3(:,i)=squeeze(DensityComp_Under1(:,ind_a(i),ind_b(i),card));                
            end
            
            %[resp, card, size(ind_a,1)]
            clear ind_a ind_b
                       
            kdata_Under2{resp, card}=tmp1;
            Traj_Under2{resp, card}=tmp2;
            DensityComp_Under2{resp, card}=tmp3;
            
            clear tmp1 tmp2 tmp3

%             tmp1=kdata_Under1(:,:,:,(resp-1)*tres+1:resp*tres,card);
%             tmp2=Traj_Under1(:,:,(resp-1)*tres+1:resp*tres,card);
%             tmp3=DensityComp_Under1(:,:,(resp-1)*tres+1:resp*tres,card);
        
%             tmp1=permute(tmp1,[1,2,4,3]);
% 
%             kdata_Under2{resp, card}=reshape(tmp1,[nx,size(kdata_Under1,2)*size(tmp1,3),nc]);
%             Traj_Under2{resp, card}=reshape(tmp2,[nx,size(kdata_Under1,2)*size(tmp2,3)]);
%             DensityComp_Under2{resp, card}=reshape(tmp3,[nx,size(kdata_Under1,2)*size(tmp3,3)]);
        end       
    end
else
    disp('not enough spokes per image for this number of cardiac phases')
end


% % headcount
% sm =0;
%   for resp = 1:nResp
%         for card = 1:nCard
%             sm = sm + size(tmp1,2);
%         end
%   end
%   


end%of function