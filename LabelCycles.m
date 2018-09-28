function [ finalLabels, para ] = LabelCycles( Cardiac_Signal, para , condition)
% labels cardiac cycle by length
  vv = zeros(size(para.ES_index));
  vh = zeros(size(para.ES_index));
  hv = zeros(size(para.ES_index));
  cycleLabels =zeros(size(para.ES_index));
  finalLabels =zeros(size(para.ES_index));
  
 % classes = [10, 15, 20, 25];  
  %thresh = 1/para.HeartFS/para.TR + 4
        
  thresh = 5;
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
  for ii=1:(size(para.ED_index,1) -3)
    vv(ii+1) = para.ES_index(ii+1)-para.ES_index(ii);    
    hv(ii+1) = para.ES_index(ii+1)-para.ED_index(ii+1);
    vh(ii+1) = para.ED_index(ii+2)-para.ES_index(ii+1);
  end
 
  if (length(para.ED_index)>para.ntres)
    vh(para.ntres) = para.ED_index(para.ntres)-para.ES_index(para.ntres -1);
  else
    vh(para.ntres) = length(Cardiac_Signal)-para.ES_index(length(para.ES_index));
  end
  
  H = histogram(vh,[1:50]);
  
  %find max position
  maxBin = find(H.Values==max(H.Values))
  %lower cutoff
  lowerThresh = floor(maxBin*0.5)
  
  binThresh = 20; %15
  counter = 0;
  
  smoothHisto = smooth(H.Values, 0.1, 'loess');
  [peaks, indpeaks] = findpeaks(smoothHisto);
  newpeaks = []
  newindpeaks = []
  for p=1:length(peaks)
      if peaks(p) > 5
          newpeaks = [newpeaks, [peaks(p)]]
          newindpeaks = [newindpeaks, [indpeaks(p)]]
      end
  end
  indpeaks = newindpeaks
  peaks = newpeaks
  
  if strcmp(condition,'pvc') % short, normal           
           cls = 2;
           if length(peaks)>1
              if peaks(1) > peaks(2)
                  secondMaxBin = indpeaks(2);
                  dist1 = abs(secondMaxBin - maxBin)*0.5;
                  upperThresh = secondMaxBin+3;
                  for ind = int8(lowerThresh):int8(maxBin+dist1)
                      counter = counter + H.Values(ind);
                      pos = find(vh==ind);
                      cycleLabels(pos) = 1;
                  end
                  for ind = int8(maxBin+dist1)+1:upperThresh
                      counter = counter + H.Values(ind);
                      pos = find(vh==ind);
                      cycleLabels(pos) = 2;
                  end
              else
                  secondMaxBin = indpeaks(1);
                  dist1 = abs(secondMaxBin - maxBin)*0.5;
                  upperThresh = maxBin+4;
                  lowerThresh = floor(secondMaxBin*0.5)
                  for ind = int8(lowerThresh):int8(secondMaxBin+dist1)
                      counter = counter + H.Values(ind);
                      pos = find(vh==ind);
                      cycleLabels(pos) = 1;
                  end
                  for ind = int8(secondMaxBin+dist1)+1:upperThresh
                      counter = counter + H.Values(ind);
                      pos = find(vh==ind);
                      cycleLabels(pos) = 2;
                  end
              end             
           end           
  else
     if strcmp(condition,'2block') % short, longer, long           
%            cls = 4;
%            for ind = int8(lowerThresh):10
%                counter = counter + H.Values(ind);
%                pos = find(vh==ind);
%                cycleLabels(pos) = 1;
%            end
%            for ind = 10:20
%                counter = counter + H.Values(ind);
%                pos = find(vh==ind);
%                cycleLabels(pos) = 2;
%            end
%            for ind = 20:30
%                counter = counter + H.Values(ind);
%                pos = find(vh==ind);
%                cycleLabels(pos) = 3;
%            end
%            for ind = 30:40
%                counter = counter + H.Values(ind);
%                pos = find(vh==ind);
%                cycleLabels(pos) = 4;
%            end
           if length(peaks)>2
              secondMaxBin = indpeaks(2);
              thirdMaxBin = indpeaks(3);
           end
           dist1 = (secondMaxBin - maxBin)*0.5;
           dist2 = (thirdMaxBin - secondMaxBin)*0.5;
           upperThresh = thirdMaxBin+3;
           for ind = int8(lowerThresh):int8(maxBin+dist1)
               counter = counter + H.Values(ind);
               pos = find(vh==ind);
               cycleLabels(pos) = 1;
           end
           for ind = int8(maxBin+dist1)+1:int8(secondMaxBin+dist2)
               counter = counter + H.Values(ind);
               pos = find(vh==ind);
               cycleLabels(pos) = 2;
           end
           for ind = int8(secondMaxBin+dist2)+1:upperThresh
               counter = counter + H.Values(ind);
               pos = find(vh==ind);
               cycleLabels(pos) = 3;
           end
     else
       if strcmp(condition,'afib') % random beats             
           %upper cutoff
           upperThresh = 30;% maxBin*1.4;         
           cls = 5;
           ind = maxBin+1;
           % put all the bins after the peak in separate classes
           while (ind <= upperThresh)
               counter = H.Values(ind);
               pos = find(vh==ind);
               cycleLabels(pos) = cls;
               while (counter < binThresh && ind <= upperThresh)
                   ind = ind +1 ;
                   counter = counter + H.Values(ind);
                   pos = find(vh==ind);
                   cycleLabels(pos) = cls;
               end
               ind = ind + 1;
               cls = cls + 1;
           end
           if (counter < binThresh)
               for ind = 1:length(cycleLabels)
                   if (cycleLabels(ind) == cls-1)
                       cycleLabels(ind) = cls-2;
                   end
               end
           end
           
           % heuristically classify the bins before the peak
           ind = maxBin;
           cls = 4;
           while (ind >= lowerThresh)
               counter = H.Values(ind);
               pos = find(vh==ind);
               cycleLabels(pos) = cls;
               while (counter < binThresh && ind >= lowerThresh)
                   ind = ind -1 ;
                   counter = counter + H.Values(ind);
                   pos = find(vh==ind);
                   cycleLabels(pos) = cls;
               end
               ind = ind - 1;
               cls = cls - 1;
           end
           if (counter < binThresh)
               for ind = 1:length(cycleLabels)
                   if (cycleLabels(ind) == cls+1)
                       cycleLabels(ind) = 0;
                   end
               end
           end   
       else
          cycleLabels(:) = 1;
       end
     end
  end
      
%  for ii=2:(para.ntres)
%     if (vv(ii) > 15)
%        if (vv(ii) < 25 )        
%          cycleLabels(ii) = 1;
%        end
%     end
%     if (vh(ii-1) > thresh)
%       if (vh(ii-1) > thresh + 1)
%          z(ii) = 2;  
%       else
%          cycleLabels(ii) = 1;
%       end
%     end
%       j = 4;
%       while (vh(ii-1) > j && j<8)
%          j = j+1;    
%       end
%       cycleLabels(ii) = j-3;
%       if (vh(ii-1) < 20)          
%         if (vh(ii-1) < 8)
%           cycleLabels(ii) = 1;
%         elseif (vh(ii-1) < 9)
%           cycleLabels(ii) = 2; 
%         elseif (vh(ii-1) < 10)
%           cycleLabels(ii) = 3; 
%         elseif (vh(ii-1) < 11)
%           cycleLabels(ii) = 4;
%         elseif (vh(ii-1) < 12)
%           cycleLabels(ii) = 5;       
%         else
%           cycleLabels(ii) = 6;
%         end
%       end

%       for j=1:length(classes)
%           if (valley_size(ii) <= classes(j)+ thresh) && (valley_size(ii) >= classes(j) - thresh)
%               if cycleLabels(ii) == 0
%                   cycleLabels(ii) = j;
%               else
%                   cycleLabels(ii) = cycleLabels(ii) + 0.5;
%               end
%           end
%       end
%  end
%  if (para.ED_index(1) > thresh)
%    cycleLabels(1) = 1;
%  end
  cycleLabels(:) = cycleLabels(:) + 1 - min(cycleLabels(find(cycleLabels>0)));
  cycleLabels(find(cycleLabels<0)) = 0;
  
  for ii=2:(length(cycleLabels))
      if (cycleLabels(ii-1) > 0)
          finalLabels(ii) = cycleLabels(ii-1);
      end
  end

  cycle_len = vh(3:length(vv));
  cycle_prev_len = vh(2:length(vv)-1);
  
  jitterammount = 0.5;
  jitterx = 2*(rand(size(cycle_len))-0.5)*jitterammount;
  jittery = 2*(rand(size(cycle_len))-0.5)*jitterammount;
  
  figure, hold on
  scatter(cycle_len+jitterx,cycle_prev_len+jittery)
  %plot(cycle_len,cycle_prev_len)
  hold off
  
  figure, H = histogram(vh,[1:50]);
  figure, H1 = histogram(cycleLabels,max(cycleLabels)+1);

  para.sys_size = hv;
  para.dia_size = vh;
  
end

