function mask_recon_GRASP = rr_tcm_createCoilMasks(kdata_Under, Traj_Under, DensityComp_Under, b1)

    param.E=MCNUFFT_uncomb(Traj_Under,DensityComp_Under,b1);
    param.y=double(squeeze(kdata_Under));
    recon_GRASP = param.E'*param.y; % goal: get coil images individually

    figure,imagescn(abs(recon_GRASP),[0 .0001],[],[],4)

   %%
med_recon_GRASP = medfilt1(abs(recon_GRASP),3,[],4);
averaged_recon_GRASP = squeeze(mean(recon_GRASP, 4)); % --> think about frequency domain
entropies = zeros(size(recon_GRASP,3),1);
entropies1 = zeros(size(recon_GRASP,3),1);


for ch = 1:size(recon_GRASP,3)
    averaged_img = mat2gray(abs(averaged_recon_GRASP(:,:,ch)), [min(min(abs(averaged_recon_GRASP(:,:,ch)))) max(max(abs(averaged_recon_GRASP(:,:,ch))))]);
    entropies(ch,1) = entropy(averaged_img);
    entropy(averaged_img)
    entropies1(ch,1) = sum(sum(entropyfilt(averaged_img)));
    sum(sum(entropyfilt(averaged_img)))
    
    averaged_img = medfilt2(averaged_img, [5 5]);
    averaged_img = imgaussfilt(averaged_img, 4);
    
    border_img = averaged_img; % maybe justs crop data
    border_size = 100;
    for x = 1:size(averaged_img,1)
        for y = 1:size(averaged_img,2)
            averaged_img(x,y) = sqrt(averaged_img(x,y));
            if (x < border_size || x > size(averaged_img,1)-border_size) ...
                    || (y < border_size || y > size(averaged_img,1)-border_size)
                border_img(x,y) = 0;
            else
                border_img(x,y) = 1;
            end
        end
    end

      % compute automatic threshold
      [counts, x] = imhist(averaged_img,128);
      thresh = otsuthresh(counts);

      % compute connected components (cc)
      bw = im2bw(averaged_img, thresh);
      bw = imfill(bw, 'holes');
      [labeledImage, numberOfBlobs] = bwlabel(bw);
     
      se = strel('octagon', 6);
      
      % get largest cc
      largestCC = extractNLargestCC(bw, 1);
      largestCC = imclose(largestCC, se);
      largestCC = largestCC.*border_img; 
            
      D = bwdist(largestCC);       
      cD = imcomplement(D)+max(D(:));
      cD = imopen(cD, se); 
      cD = medfilt2(cD);
      mask_recon_GRASP(:,:,ch) = mat2gray(cD, [double(floor(max(cD(:))*.3)),double(floor(max(cD(:))*.95))]);
      
      centerGaussian = fspecial('gaussian', [size(bw,1) size(bw,2)], size(bw,1)/4);
      centerGraussianImg = mat2gray(centerGaussian, [0, max(max(centerGaussian))*3/4]);
      
      mask_recon_GRASP(:,:,ch) = max(centerGraussianImg, mask_recon_GRASP(:,:,ch));
      mask_recon_GRASP(:,:,ch) = medfilt2(mask_recon_GRASP(:,:,ch));

end

    figure,imagescn(squeeze(mask_recon_GRASP),[0 1],[],[],4)
end

