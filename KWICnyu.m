% calculates KWIC mask
% firstRingSize can be any arbitrary number

function [mask, dcf] = KWICnyu(Nread, Nproj, nyqFactor, firstRingSize, centerP, level)

if exist('nyqFactor')~=1
    nyqFactor = 1.1; % the higher, the more data is used for each time frame
end

if exist('firstRingSize')~=1
    firstRingSize = 34; % or 8/13/21/34 (fibonacci numbers guarantee near-optimal k-space sampling)
end

os = 2;
k  = 1;
KWICrad(1) = 0;

while fibonacci(k) < 2*Nproj
    KWICrad     (k+1) = fibonacci(k+1) / pi / nyqFactor * os;
    KWICprojN   (k) = fibonacci(k+1);
    k = k+1;
end

KWICprojN = [KWICprojN,firstRingSize]; % insert firstRingSize in array
KWICprojN = sort(KWICprojN,'ascend'); % get the right order again
KWICprojN(KWICprojN<firstRingSize) = [];

KWICrad(1:length(KWICrad) - length(KWICprojN)) = [] % adapt KWICrad

if exist('centerP')~=1
    centerP = floor(KWICprojN(end)/2)+1;
end

if exist('level')~=1
    level = 0.5;
end

ival_old = centerP;
maxProj = Nproj;
mask = ones(Nread,maxProj);

for ring=1:length(KWICrad),
    
    ival = -floor(KWICprojN(ring)/2):-floor(KWICprojN(ring)/2)+KWICprojN(ring)-1;
    ival = centerP+ival;
    
    overflowLeft  =   max(0,1-ival(1)); % makes sure ival starts at least at 1
    overflowRight =   max(0,ival(end)-maxProj ); % makes sure ival ends at Nproj
    
    ival =  ival + overflowLeft - overflowRight;
    
    ival = max(ival,1);
    ival = min(ival,Nproj);
    
    ival_new = ival(ival<ival_old(1) | ival>ival_old(end))
    ival_old = ival;
    
    if ring > 1
        mask(max(Nread/2-floor(KWICrad(ring-1)),1) : min(Nread/2+floor(KWICrad(ring-1)),Nread),ival_new) = repmat(1-tukeywin(size(max(Nread/2-floor(KWICrad(ring-1)),1) : min(Nread/2+floor(KWICrad(ring-1)),Nread),2),level),[1 size(ival_new,2)]);
    end
end

dcf = repmat(maxProj./sum(mask,2),[1, maxProj]);
