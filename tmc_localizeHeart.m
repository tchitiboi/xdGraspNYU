function mask = tmc_localizeHeart(imgs, HF_Index)

permuted_imgs = permute(imgs, [3 2 1]);
fft_imgs = fftshift(fft(permuted_imgs), 1);
img_fft_imgs = mat2gray(real(fft_imgs));
selected_frequency = squeeze(img_fft_imgs(HF_Index,:,:));

accum = zeros(size(imgs, 2), size(imgs, 1));

for t=1:max(size(HF_Index))-1
    accum = accum + squeeze((selected_frequency(t,:,:)-selected_frequency(t+1,:,:)).*(selected_frequency(t,:,:)-selected_frequency(t+1,:,:)));
end

map = accum;
map = medfilt2(map);
im_map = mat2gray(map, [0 max(max(map))]);

[counts, x] = imhist(im_map,128);
thresh = otsuthresh(counts);

% threshold and dilate
bw = im2bw(im_map, thresh);
se = strel('octagon',9);
bw_dilated = imdilate(bw,se);
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
closed_largestCC = imclose(largestCC,se);
closed_largestCC = imdilate(largestCC,se);
largestCC = largestCC.*border_img; 
closed_largestCC = imfill(closed_largestCC, 'holes');

%figure,imagescn(abs(ipermute(img_fft_imgs, [3 2 1])),[0 0.5],[],[],3)
figure,imagescn(abs(ipermute(closed_largestCC, [2 1])),[0 1],[],[],2)

mask = abs(ipermute(closed_largestCC, [2 1]));

end