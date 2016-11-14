function [Res_Signal_Bins, Res_Signal_P] = getRespBins(Res_Signal, nbins)
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
    while H.Values(minBin)<3
      minBin = minBin+1;  
    end
    minRes = absmin + (minBin-1)*(absmax-absmin)/100;
    
    maxBin = 100;
    while H.Values(maxBin)<3
      maxBin = maxBin-1;  
    end
    maxRes = absmin + maxBin*(absmax-absmin)/100;
    
    indx = 0;    
    for r = 1:size(Res_Signal,1)
        if (Res_Signal(r) >= minRes && Res_Signal(r) <= maxRes)
          indx = indx + 1; 
          Res_Signal_no_outliers(indx) = r;
        end        
    end
        
%     [idx,c]=kmeans(Res_Signal(Res_Signal_no_outliers),nbins);
%     index = 1:(nbins);
%     m = [c,index.'];
%     m = sortrows(m,1);
    
    guess = zeros(size(Res_Signal(Res_Signal_no_outliers))); 
    for i=1:length(guess)
        guess(i) = floor((Res_Signal(Res_Signal_no_outliers(i))-minRes)*nbins/(maxRes+0.00001-minRes)+1);
    end
    
    for i=1:nbins
        step = 0; 
        while length((find(guess==i)))<3
            step = step + 0.001;
            for j=1:length(guess)
                bin = floor((Res_Signal(Res_Signal_no_outliers(j))-minRes)*nbins/(maxRes+0.00001-minRes)+1);
                if ((Res_Signal(Res_Signal_no_outliers(j))-minRes) < i*(maxRes-minRes)/nbins+step &&...
                        (Res_Signal(Res_Signal_no_outliers(j))-minRes) > (i-1)*(maxRes-minRes)/nbins-step)
                  guess(j) = i;
                end
            end
        end
    end
    
    
    gmfit = fitgmdist(Res_Signal(Res_Signal_no_outliers), nbins, 'CovarianceType', ...
        'diagonal', 'SharedCovariance', false, 'Replicates', 1, 'Start', guess);
    idx = cluster(gmfit,Res_Signal(Res_Signal_no_outliers));
    P = posterior(gmfit,Res_Signal(Res_Signal_no_outliers));
    c = gmfit.mu;
    index = 1:(nbins);
    m = [c,index.'];
    m = sortrows(m,1);
    P1 = P;
    
    for i = 1:nbins;
        P1(:,i) = P(:,find(m(:,2)==i));
    end     
    
    Res_Signal_Bins = zeros(size(Res_Signal));
    Res_Signal_P = zeros([size(Res_Signal,1),nbins]);
    
    for r = 1:size(idx,1)
        i = idx(r);
        p = P1(r);
        Res_Signal_Bins(Res_Signal_no_outliers(r)) = find(m(:,2)==i);
        Res_Signal_P(Res_Signal_no_outliers(r),:)=P1(r,:);
    end
    
    Res_Signal_Bins = medfilt1(Res_Signal_Bins,7);
 
    
%     kmHisto = zeros(nbins,1);
%     for bin = 1:nbins
%        l = find(Res_Signal_Bins==bin);
%        kmHisto(bin)=length(l);
%     end

end

