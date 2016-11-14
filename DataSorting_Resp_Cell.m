function [kdata_Under2,Traj_Under2,DensityComp_Under2,p2]=DataSorting_Resp_Cell(kdata,Traj,DensityComp,Res_Signal,Res_Signal_P, nline,para, nResp, nCard);
%Data sorting into two dynamic dimensions, one cardiac and one respiratory

[nx,ntviews,nc]=size(kdata);
kdata_Under=zeros(nx,para.CardiacPhase*nline,nc,para.ntres);
Traj_Under=zeros(nx,para.CardiacPhase*nline,para.ntres);
DensityComp_Under=zeros(nx,para.CardiacPhase*nline,para.ntres);
Res_Signal_Under=zeros(para.CardiacPhase*nline, para.ntres);
Res_Signal_P_Under=zeros(para.CardiacPhase*nline, para.ntres, size(Res_Signal_P,2));

Res_Signal_Long = interp1( linspace(0,1,length(Res_Signal)), Res_Signal, linspace(0,1,ntviews), 'nearest');
Res_Signal_P_Long = interp2(Res_Signal_P.', (1+(1:ntviews)*(length(Res_Signal)-1)/ntviews).', (1:nResp),'nearest');
Res_Signal_P_Long = Res_Signal_P_Long.';


%Separate cardiac cycles
for ii=1:para.ntres
    kdata_Under(:,:,:,ii)=kdata(:,(para.ES_index(ii)-para.Sys)*nline+1:(para.ES_index(ii)+para.Dia)*nline,:);
    Traj_Under(:,:,ii)=Traj(:,(para.ES_index(ii)-para.Sys)*nline+1:(para.ES_index(ii)+para.Dia)*nline);
    DensityComp_Under(:,:,ii)=DensityComp(:,(para.ES_index(ii)-para.Sys)*nline+1:(para.ES_index(ii)+para.Dia)*nline);
    Res_Signal_Under(:,ii)=Res_Signal_Long((para.ES_index(ii)-para.Sys)*nline+1:(para.ES_index(ii)+para.Dia)*nline);
    Res_Signal_P_Under(:,ii,:)=Res_Signal_P_Long((para.ES_index(ii)-para.Sys)*nline+1:(para.ES_index(ii)+para.Dia)*nline,:);
end

%collapse cardiac phases
tcar=para.CardiacPhase/nCard;

kdata_Under1=zeros(nx,tcar*nline,nc,para.ntres,nCard);
Traj_Under1=zeros(nx,tcar*nline,para.ntres,nCard);
DensityComp_Under1=zeros(nx,tcar*nline,para.ntres,nCard);
Res_Signal_Under1=zeros(tcar*nline, para.ntres, nCard);
Res_Signal_P_Under1=zeros(tcar*nline, para.ntres, nCard, size(Res_Signal_P,2));

if (tcar>=1)
  for ii=1:nCard
    kdata_Under1(:,:,:,:,ii) = kdata_Under(:,(ii-1)*nline*tcar+1:ii*nline*tcar,:,:);
    Traj_Under1(:,:,:,ii) = Traj_Under(:,(ii-1)*nline*tcar+1:ii*nline*tcar,:);
    DensityComp_Under1(:,:,:,ii) = DensityComp_Under(:,(ii-1)*nline*tcar+1:ii*nline*tcar,:);
    Res_Signal_Under1(:,:,ii) = Res_Signal_Under((ii-1)*nline*tcar+1:ii*nline*tcar,:);
    Res_Signal_P_Under1(:,:,ii,:) = Res_Signal_P_Under((ii-1)*nline*tcar+1:ii*nline*tcar,:,:);
  end
else
    disp('not enough spokes per image for this number of cardiac phases')
end

%initializa cell array
kdata_Under2 = cell(nResp,nCard);
Traj_Under2 = cell(nResp,nCard);
DensityComp_Under2 = cell(nResp,nCard);
p2 = cell(nResp,nCard);
minSpokes = 50;

tres=para.ntres/nResp;
if (tres>=1)  
    for resp = 1:nResp
        for card = 1:nCard
            r = resp;
            [ind_a, ind_b] = find(Res_Signal_Under1(:,:,card) == r);
            ind_a_extra = [];
            ind_b_extra = [];
            
            if (length(ind_a)<minSpokes)   
                missing = minSpokes - length(ind_a);               
                [ind_a3, ind_b3] = find(Res_Signal_Under1(:,:,card) ~= r);
                p = zeros(length(ind_a3),2);
                for i = 1:length(ind_a3)                    
                    p(i,1) = Res_Signal_P_Under1(ind_a3(i),ind_b3(i),card,r);
                    p(i,2) = i;
                end
                p = sortrows(p,-1);
                ind_a_extra = ind_a3(p(1:min(missing,size(p,1)),2));
                ind_b_extra = ind_b3(p(1:min(missing,size(p,1)),2));
                p_extra = p(1:min(missing,size(p,1)),1);
            end 
            
            tmp1 = zeros(nx, length(ind_a)+length(ind_a_extra), nc);
            tmp2 = zeros(nx, length(ind_a)+length(ind_a_extra));
            tmp3 = zeros(nx, length(ind_a)+length(ind_a_extra));
            tmp4 = zeros(length(ind_a)+length(ind_a_extra),1);
                
            for i = 1:length(ind_a)
                tmp1(:,i,:)=squeeze(kdata_Under1(:,ind_a(i),:,ind_b(i),card));
                tmp2(:,i)=squeeze(Traj_Under1(:,ind_a(i),ind_b(i),card));
                tmp3(:,i)=squeeze(DensityComp_Under1(:,ind_a(i),ind_b(i),card));     
                tmp4(i)=squeeze(Res_Signal_P_Under1(ind_a(i),ind_b(i),card,r));
            end
            
            for i = 1:length(ind_a_extra)
                tmp1(:,i+length(ind_a),:)=squeeze(kdata_Under1(:,ind_a_extra(i),:,ind_b_extra(i),card));
                tmp2(:,i+length(ind_a))=squeeze(Traj_Under1(:,ind_a_extra(i),ind_b_extra(i),card));
                tmp3(:,i+length(ind_a))=squeeze(DensityComp_Under1(:,ind_a_extra(i),ind_b_extra(i),card)); 
                tmp4(i+length(ind_a))=p_extra(i);
            end
            
            clear ind_a ind_b
                       
            kdata_Under2{resp, card} = tmp1;
            Traj_Under2{resp, card} = tmp2;
            DensityComp_Under2{resp, card} = tmp3;
            p2{resp, card} = tmp4;
            
            clear tmp1 tmp2 tmp3 tmp4

        end       
    end
else
    disp('not enough spokes per image for this number of cardiac phases')
end

end%of function