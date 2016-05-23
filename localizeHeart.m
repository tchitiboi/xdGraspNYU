%Rebecca Ramb
%2016-04-08 : more precise localization of the heart for improved
%derivation of cardiac cycle
%2016-05-18 : 

function out_dilated = localizeHeart(data,threshold)

if nargin < 2
    threshold(1) = 0.2;
    threshold(2) = 0.2;
end
out_tmp = complex(zeros(size(data)));
nt = size(out_tmp,3);
for tt=1:nt
tmp = abs(data(:,:,tt));
tmp = tmp./max(tmp(:));
BW = tmp;
BW(tmp<threshold(1))=0;
BW(tmp>=threshold(1))=1;
figure,imagescn(abs(BW),[0 1],[],[],3)
CC = bwconncomp(BW);
numPixels = cellfun(@numel,CC.PixelIdxList);
[~,idx] = max(numPixels);
BW(CC.PixelIdxList{idx}) = 1000;
BW(BW<1000)=0;
BW(BW==1000)=1;
out_tmp(:,:,tt) = BW;
end
out = mean(out_tmp,3).^2;
out(out<threshold(2))=0;
CC = bwconncomp(out);
numPixels = cellfun(@numel,CC.PixelIdxList);
[~,idx] = max(numPixels);
out(CC.PixelIdxList{idx}) = 1000;
out(out<1000)=0;
out(out==1000)=1;
se = offsetstrel('ball',25,25);
out_dilated = imdilate(out,se); out_dilated = out_dilated-min(out_dilated(:));
out_dilated = (out_dilated)./max(out_dilated(:));
out_dilated(out_dilated>0)=1;
figure, imagesc(out_dilated)

end