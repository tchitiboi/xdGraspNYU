function mask_recon_GRASP = createCoilMasks(kdata_Under, Traj_Under, DensityComp_Under, b1)

    param.E=MCNUFFT_GetCoils(Traj_Under,DensityComp_Under,b1);
    param.y=double(squeeze(kdata_Under));
    recon_GRASP = param.E'*param.y;

    figure,imagescn(abs(recon_GRASP),[0 .0001],[],[],4)

    med_recon_GRASP = medfilt1(real(recon_GRASP),3,[],4);
    averaged_recon_GRASP = squeeze(mean(recon_GRASP, 4));
    entropies = zeros(size(recon_GRASP,3),1);
    entropies1 = zeros(size(recon_GRASP,3),1);

    for ch = 1:size(recon_GRASP,3)  
      averaged_img = mat2gray(abs(averaged_recon_GRASP(:,:,ch)), [min(min(abs(averaged_recon_GRASP(:,:,ch)))) max(max(abs(averaged_recon_GRASP(:,:,ch))))]);
      entropies(ch,1) = entropy(averaged_img);
      entropy(averaged_img)
      entropies1(ch,1) = sum(sum(entropyfilt(averaged_img)));
      sum(sum(entropyfilt(averaged_img)))

      if entropies(ch,1) < 5
          averaged_img = medfilt2(averaged_img, [5 5]);
          averaged_img = imgaussfilt(averaged_img, 4);

          border_img = averaged_img;
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
          %thresh = otsuthresh(counts);
 thresh = graythresh(counts);
          % compute connected components (cc)
          bw = im2bw(averaged_img, thresh);
          bw = imfill(bw, 'holes');
          [labeledImage, numberOfBlobs] = bwlabel(bw);
          %blobMeasurements = regionprops(labeledImage, 'area', 'Centroid');

          % get largest cc
          largestCC = extractNLargestCC(bw, 1);
          largestCC = largestCC.*border_img; 
          se = strel('octagon', 6);
          D = bwdist(largestCC);       
          cD = imcomplement(D);
          cD = imopen(cD, se);      

          mask_recon_GRASP(:,:,ch) = real(mat2gray(cD, [-300,-50]));

      else
          mask_recon_GRASP(:,:,ch) = real(zeros(size(recon_GRASP,1), size(recon_GRASP,1)));
      end
    end

    mask_recon_GRASP = mask_recon_GRASP.^2;
    for ch = 1:size(recon_GRASP,3) 
        for y = 1:size(recon_GRASP,2)
            for x = 1:size(recon_GRASP,1)
               if mask_recon_GRASP(x,y,ch)<0.01 
                  mask_recon_GRASP(x,y,ch)= 0.01;
               end
            end
        end
    end

    figure,imagescn(squeeze(mask_recon_GRASP),[0 1],[],[],4)
end

