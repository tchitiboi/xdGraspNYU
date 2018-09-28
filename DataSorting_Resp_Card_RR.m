function [kdata_Under2,Traj_Under2,DensityComp_Under2,p2]=DataSorting_Resp_Card_RR(kdata,Traj,DensityComp,Res_Signal,Res_Signal_P, cycleLabels, labels, nline,para, nResp, nCard);
%Data sorting into two dynamic dimensions, one cardiac and one respiratory
%para.ntres = length(find(cycleLabels==1));

[nx,ntviews,nc]=size(kdata);
% kdata_Under=zeros(nx,para.CardiacPhase*nline,nc,para.ntres);
% Traj_Under=zeros(nx,para.CardiacPhase*nline,para.ntres);
% DensityComp_Under=zeros(nx,para.CardiacPhase*nline,para.ntres);
% Res_Signal_Under=zeros(para.CardiacPhase*nline, para.ntres);
% Res_Signal_P_Under=zeros(para.CardiacPhase*nline, para.ntres, size(Res_Signal_P,2));

Res_Signal_Long = interp1( linspace(0,1,length(Res_Signal)), Res_Signal, linspace(0,1,ntviews), 'nearest');
Res_Signal_P_Long = interp2(Res_Signal_P.', (1+(1:ntviews)*(length(Res_Signal)-1)/ntviews).', (1:nResp),'nearest');
Res_Signal_P_Long = Res_Signal_P_Long.';


%Separate cardiac cycles
% for ii=1:para.ntres
%     if abs(cycleLabels(ii)-label) < 1
%         kdata_Under(:,:,:,ii)=kdata(:,(para.ES_index(ii)-para.Sys)*nline+1:(para.ES_index(ii)+para.Dia)*nline,:);
%         Traj_Under(:,:,ii)=Traj(:,(para.ES_index(ii)-para.Sys)*nline+1:(para.ES_index(ii)+para.Dia)*nline);
%         DensityComp_Under(:,:,ii)=DensityComp(:,(para.ES_index(ii)-para.Sys)*nline+1:(para.ES_index(ii)+para.Dia)*nline);
%         Res_Signal_Under(:,ii)=Res_Signal_Long((para.ES_index(ii)-para.Sys)*nline+1:(para.ES_index(ii)+para.Dia)*nline);
%         Res_Signal_P_Under(:,ii,:)=Res_Signal_P_Long((para.ES_index(ii)-para.Sys)*nline+1:(para.ES_index(ii)+para.Dia)*nline,:);
%     end
% end

kdata_cycle = cell(para.ntres,labels);
Traj_cycle = cell(para.ntres,labels);
DensityComp_cycle = cell(para.ntres,labels);
Res_Signal_cycle = cell(para.ntres,labels);
Res_Signal_P_cycle = cell(para.ntres,labels);

for ii=1:para.ntres-1
  for label=1:labels
    if abs(cycleLabels(ii)-label) < 1
        kdata_cycle{ii,label}=kdata(:,(para.ES_index(ii)-para.sys_size(ii))*nline+1:(para.ES_index(ii)+para.dia_size(ii))*nline,:);
        Traj_cycle{ii,label}=Traj(:,(para.ES_index(ii)-para.sys_size(ii))*nline+1:(para.ES_index(ii)+para.dia_size(ii))*nline);
        DensityComp_cycle{ii,label}=DensityComp(:,(para.ES_index(ii)-para.sys_size(ii))*nline+1:(para.ES_index(ii)+para.dia_size(ii))*nline);
        Res_Signal_cycle{ii,label}=Res_Signal_Long((para.ES_index(ii)-para.sys_size(ii))*nline+1:(para.ES_index(ii)+para.dia_size(ii))*nline);
        Res_Signal_P_cycle{ii,label}=Res_Signal_P_Long((para.ES_index(ii)-para.sys_size(ii))*nline+1:(para.ES_index(ii)+para.dia_size(ii))*nline,:);  
    end
  end
end

%collapse cardiac phases
tcar=para.CardiacPhase/nCard;
% 
% kdata_Under1=zeros(nx,floor(tcar*nline),nc,para.ntres,nCard);
% Traj_Under1=zeros(nx,floor(tcar*nline),para.ntres,nCard);
% DensityComp_Under1=zeros(nx,floor(tcar*nline),para.ntres,nCard);
% Res_Signal_Under1=zeros(floor(tcar*nline), para.ntres, nCard);
% Res_Signal_P_Under1=zeros(floor(tcar*nline), para.ntres, nCard, size(Res_Signal_P,2));

kdata_percycle = cell(para.ntres,nCard,labels);
Traj_percycle = cell(para.ntres,nCard,labels);
DensityComp_percycle = cell(para.ntres,nCard,labels);
Res_Signal_percycle = cell(para.ntres,nCard,labels);
Res_Signal_P_percycle = cell(para.ntres,nCard,labels);

%binsize = floor(nline*tcar);

if (tcar>=0.5)
  for ii=1:nCard
    for j=1:para.ntres
      for label=1:labels
          if ~isempty(kdata_cycle{j,label})
            %j=j
            binsize = floor((para.sys_size(j) + para.dia_size(j))/nCard*nline);
            maxLen = size(kdata_cycle{j,label},2);
            kdata_percycle{j,ii,label} = kdata_cycle{j,label}(:,(ii-1)*binsize+1:min(ii*binsize,maxLen),:);
            Traj_percycle{j,ii,label} = Traj_cycle{j,label}(:,(ii-1)*binsize+1:min(ii*binsize,maxLen));
            DensityComp_percycle{j,ii,label} = DensityComp_cycle{j,label}(:,(ii-1)*binsize+1:min(ii*binsize,maxLen));
            Res_Signal_percycle{j,ii,label} = Res_Signal_cycle{j,label}((ii-1)*binsize+1:min(ii*binsize,maxLen));
            Res_Signal_P_percycle{j,ii,label} = Res_Signal_P_cycle{j,label}((ii-1)*binsize+1:min(ii*binsize,maxLen),:);
          end  
      end
    end
  end
else
    disp('not enough spokes per image for this number of cardiac phases')
end

%initializa cell array
kdata_Under2 = cell(nResp,nCard,labels);
Traj_Under2 = cell(nResp,nCard,labels);
DensityComp_Under2 = cell(nResp,nCard,labels);
p2 = cell(nResp,nCard,labels);
minSpokes = 25;

tres=para.ntres/nResp;
if (tres>=1)  
  for label=1:labels
    for resp = 1:nResp
        for card = 1:nCard
            r = resp;
            ind_a=[];
            ind_b=[];
            for j = 1:para.ntres
               ind_b = [ind_b,find(Res_Signal_percycle{j,card,label} == r)];
               ind_a = [ind_a,repmat(j,1,length(find(Res_Signal_percycle{j,card,label} == r)))];
            end
            ind_a_extra = [];
            ind_b_extra = [];
            
            if (length(ind_a)<minSpokes)   
                missing = minSpokes - length(ind_a);               
                %[ind_a3, ind_b3] = find(Res_Signal_percycle{j,card}(:,:) ~= r);
                ind_a3=[];
                ind_b3=[];
                for j = 1:para.ntres
                  ind_b3 = [ind_b3,find(Res_Signal_percycle{j,card,label} ~= r)];
                  ind_a3 = [ind_a3,repmat(j,1,length(find(Res_Signal_percycle{j,card,label} ~= r)))];
                end
                p = zeros(length(ind_a3),2);
                for i = 1:length(ind_a3)                    
                    p(i,1) = Res_Signal_P_percycle{ind_a3(i),card,label}(ind_b3(i),r);
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
                tmp1(:,i,:)=squeeze(kdata_percycle{ind_a(i),card,label}(:,ind_b(i),:));
                tmp2(:,i)=squeeze(Traj_percycle{ind_a(i),card,label}(:,ind_b(i)));
                tmp3(:,i)=squeeze(DensityComp_percycle{ind_a(i),card,label}(:,ind_b(i)));     
                tmp4(i)=squeeze(Res_Signal_P_percycle{ind_a(i),card,label}(ind_b(i),r));
            end
            
            for i = 1:length(ind_a_extra)
                tmp1(:,i+length(ind_a),:)=squeeze(kdata_percycle{ind_a_extra(i),card,label}(:,ind_b_extra(i),:));
                tmp2(:,i+length(ind_a))=squeeze(Traj_percycle{ind_a_extra(i),card,label}(:,ind_b_extra(i)));
                tmp3(:,i+length(ind_a))=squeeze(DensityComp_percycle{ind_a_extra(i),card,label}(:,ind_b_extra(i))); 
                tmp4(i+length(ind_a))=p_extra(i);
            end
            
            clear ind_a ind_b
                       
            kdata_Under2{resp,card,label} = tmp1;
            Traj_Under2{resp,card,label} = tmp2;
            DensityComp_Under2{resp,card,label} = tmp3;
            p2{resp,card,label} = tmp4;
            
            clear tmp1 tmp2 tmp3 tmp4

        end       
    end
  end
else
    disp('not enough spokes per image for this number of resp phases')
end

end%of function