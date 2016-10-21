function Res_Signal_Bins = getRespBins(Res_Signal, nbins)
%GETRESPBINS Summary of this function goes here
%   Detailed explanation goes here

    m = mean(Res_Signal,1);
    sd = std(Res_Signal,1);
    %minRes = m - 1.96*sd;
    absmin = min(Res_Signal);
    absmax = max(Res_Signal);
    %maxRes = m + 1.96*sd;
    
    H = histogram(Res_Signal,100);
    minBin = 1;
    while H.Values(minBin)<10
      minBin = minBin+1;  
    end
    minRes = absmin + (minBin-1)*(absmax-absmin)/100;
    
    maxBin = 100;
    while H.Values(maxBin)<10
      maxBin = maxBin-1;  
    end
    maxRes = absmin + maxBin*(absmax-absmin)/100;
    
    Res_Signal_no_outliers = zeros(size(Res_Signal));
    
    for r = 1:size(Res_Signal,1)
        if (Res_Signal(r) < minRes || Res_Signal(r) > maxRes)
          Res_Signal_no_outliers(r) = -100;
        else
          Res_Signal_no_outliers(r) = Res_Signal(r);
        end        
    end
    
    [idx,c]=kmeans(Res_Signal_no_outliers,nbins+1);
    index = 1:(nbins+1);
    m = [c,index.'];
    m = sortrows(m,1);
    
    for r = 1:size(Res_Signal,1)
        i = idx(r);
        Res_Signal_Bins(r) = find(m(:,2)==i)-1;
    end
    
    Res_Signal_Bins = medfilt1(Res_Signal_Bins,7);
    
%     binsz = (maxRes-minRes)/nbins;
%     Res_Signal_Bins = zeros(size(Res_Signal));
% 
%     for r = 1:size(Res_Signal,1)
%         if (Res_Signal(r) < minRes || Res_Signal(r) > maxRes)
%              Res_Signal_Bins(r) = 0;
%         else
%              Res_Signal_Bins(r) = 1 + floor((Res_Signal(r)-minRes)/binsz); 
%         end
%     end
%     
%     for bin = 0:nbins
%        l = find(Res_Signal_Bins==bin);
%        [bin,length(l)]
%     end

end

