function [ cycleLabels, para ] = LabelCycles( Cardiac_Signal, para )
% labels cardiac cycle by length
  vv = zeros(size(para.ES_index));
  vh = zeros(size(para.ES_index));
  hv = zeros(size(para.ES_index));
  cycleLabels =zeros(size(para.ES_index));
  
 % classes = [10, 15, 20, 25];
  
  %thresh = 1/para.HeartFS/para.TR + 4
  thresh = 8;
  ii = 1;
  while ii < min(size(para.ED_index,1),size(para.ES_index,1))
    if (para.ED_index(ii+1)-para.ES_index(ii)<0)
      para.ED_index(ii+1) = []
    end
    ii = ii+1
  end
  
  para.ntres = length(para.ES_index);
  
  vv(1) = para.ES_index(1);
  hv(1) = para.ES_index(1)-para.ED_index(1);
  vh(1) = para.ED_index(2)-para.ES_index(1);
  for ii=1:(size(para.ED_index,1) -2)
    vv(ii+1) = para.ES_index(ii+1)-para.ES_index(ii);    
    hv(ii+1) = para.ES_index(ii+1)-para.ED_index(ii+1);
    vh(ii+1) = para.ED_index(ii+2)-para.ES_index(ii+1);
  end
  
  vv(size(para.ED_index,1)) = para.ES_index(size(para.ES_index,1))-para.ES_index(size(para.ED_index,1)-1);    
  hv(size(para.ED_index,1)) = para.ES_index(size(para.ES_index,1))-para.ED_index(size(para.ED_index,1));
  
  if (length(para.ED_index)>para.ntres)
    vh(para.ntres) = para.ED_index(para.ntres)-para.ES_index(para.ntres -1);
  else
    vh(para.ntres) = length(Cardiac_Signal)-para.ES_index(length(para.ES_index));
  end
  
  for ii=2:(para.ntres)
%     if (vv(ii) > 15)
%        if (vv(ii) < 25 )        
%          cycleLabels(ii) = 1;
%        end
%     end
%     if (vh(ii-1) > thresh)
%       if (vh(ii-1) > thresh + 1)
%          cycleLabels(ii) = 2;  
%       else
%          cycleLabels(ii) = 1;
%       end
%     end
%       j = 4;
%       while (vh(ii-1) > j && j<8)
%          j = j+1;    
%       end
%       cycleLabels(ii) = j-3;
      if (vh(ii-1) < 20)          
        if (vh(ii-1) < 8)
          cycleLabels(ii) = 1;
        elseif (vh(ii-1) < 9)
          cycleLabels(ii) = 2; 
        elseif (vh(ii-1) < 10)
          cycleLabels(ii) = 3; 
        elseif (vh(ii-1) < 11)
          cycleLabels(ii) = 4;
        elseif (vh(ii-1) < 12)
          cycleLabels(ii) = 5;       
        else
          cycleLabels(ii) = 6;
        end
      end

%       for j=1:length(classes)
%           if (valley_size(ii) <= classes(j)+ thresh) && (valley_size(ii) >= classes(j) - thresh)
%               if cycleLabels(ii) == 0
%                   cycleLabels(ii) = j;
%               else
%                   cycleLabels(ii) = cycleLabels(ii) + 0.5;
%               end
%           end
%       end
  end
  if (para.ED_index(1) > thresh)
    cycleLabels(1) = 1;
  end
  
  cycle_len = vv(2:length(vv));
  cycle_prev_len = vv(1:length(vv)-1);
  
  figure, hold on
  scatter(cycle_len,cycle_prev_len)
  plot(cycle_len,cycle_prev_len)
  hold off
  
  figure, H = histogram(vh,[1:50]);
  figure, H1 = histogram(cycleLabels,6);

  para.sys_size = hv;
  para.dia_size = vh;
  
end

