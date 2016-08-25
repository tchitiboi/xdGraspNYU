function mask = tmc_localizeHeart(imgs, HF_Index)

for t = 1:size(imgs,3)
  smooth_imgs(:,:,t) = imgaussfilt(imgs(:,:,t),3);
end

%smooth_imgs = medfilt1(smooth_imgs, 3, [], 3);

permuted_imgs = permute(smooth_imgs, [3 2 1]);
fft_imgs = fftshift(fft(permuted_imgs), 1);

for x = 1:size(smooth_imgs,1)
    for y = 1:size(smooth_imgs,2)
        fft_imgs(:,y,x) = fft_imgs(:,y,x)/max(fft_imgs(:,y,x));
    end
end

img_fft_imgs = mat2gray(real(fft_imgs));

selected_frequency = squeeze(img_fft_imgs(HF_Index,:,:));
accum = zeros(size(imgs, 2), size(imgs, 1));

for t=1:max(size(HF_Index))-1
    accum = accum + (medfilt2(squeeze(selected_frequency(t,:,:)))-medfilt2(squeeze(selected_frequency(t+1,:,:)))).^2;
end


map = accum;
map = medfilt2(map);
im_map = mat2gray(map, [0 max(max(map))]);

[counts, x] = imhist(im_map,256);
thresh = otsuthresh(counts)
thresh = thresh;% * 2/3.0
%thresh = graythresh(counts)

% threshold and dilate
bw = im2bw(im_map, thresh);
se = strel('octagon',6);
bw_dilated = imdilate(bw,se);
se = strel('octagon',6);
bw_dilated = imclose(bw_dilated,se);
bw_dilated = imfill(bw_dilated, 'holes');

% get largest cc
[labeledImage, numberOfBlobs] = bwlabel(bw_dilated);
largestCC = extractNLargestCC(bw_dilated, 1);

border_img = largestCC;
border_size = 30;
for x = 1:size(largestCC,1)
  for y = 1:size(largestCC,2)
      if (x < border_size || x > size(largestCC,1)-border_size) ...
         || (y < border_size || y > size(largestCC,1)-border_size)
          border_img(x,y) = 0;
      else
          border_img(x,y) = 1;  
      end
  end
end

largestCC = largestCC.*border_img; 
se = strel('octagon',9);
closed_largestCC = imclose(largestCC,se);
closed_largestCC = imdilate(largestCC,se);
closed_largestCC = imclose(largestCC,se);
%largestCC = largestCC.*border_img; 
closed_largestCC = imfill(closed_largestCC, 'holes');

%figure,imagescn(abs(ipermute(img_fft_imgs, [3 2 1])),[0 0.5],[],[],3)
figure,imagescn(abs(ipermute(closed_largestCC, [2 1])),[0 1],[],[],2)
figure,imagescn(abs(accum),[0 1],[],[],2)

mask = abs(ipermute(closed_largestCC, [2 1]));

end