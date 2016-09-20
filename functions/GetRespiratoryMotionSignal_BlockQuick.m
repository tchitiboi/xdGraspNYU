function [Res_Signal1,para]=GetRespiratoryMotionSignal_BlockQuick(para,maskHeart,ResSort,recon_Res);

%dpad = floor((size(recon_Res,1) - size(maskHeart,1))/2);
%mask = padarray(maskHeart,[dpad, dpad],0,'both');
%if size(recon_Res,1) - size(mask,1) > 0
%  mask = padarray(mask,[1, 1],0,'pre');
%end
% se = strel('octagon',3);
% mask = imdilate(maskHeart,se);

recon_Res = medfilt1(abs(recon_Res),3,[],3);

for t = 1:size(recon_Res,3)
  recon_Res(:,:,t) = imgaussfilt(abs(recon_Res(:,:,t)),2);
end

[nx,ny,nt]=size(recon_Res);
%recon_Res = recon_Res .* repmat(imcomplement(maskHeart),[1 1 size(recon_Res,3)]);

border_img = squeeze(recon_Res(:,:,1));
border_size = 35;
for x = 1:size(recon_Res,1)
  for y = 1:size(recon_Res,2)
      if (x < border_size || x > size(recon_Res,1)-border_size) ...
         || (y < border_size || y > size(recon_Res,2)-border_size)
          border_img(x,y) = 1;
      else
          border_img(x,y) = 0;  
      end
  end
end

recon_Res = recon_Res .* repmat(border_img,[1 1 size(recon_Res,3)]);

TR=para.TR*2;
time = TR:TR:nt*TR;
F_S = 1/TR;F_X = 0:F_S/(nt-1):F_S;
F_X=F_X-F_S/2;  %%% frequency after FFT of the motion signal
if mod(nt,2)==0
    F_X=F_X+F_X(nt/2);
end
FR_Index=find(F_X<para.HF_R & F_X>para.LF_R);
FC_Index=find(F_X<para.HF_H & F_X>para.LF_H);

[nx,ny,nt]=size(recon_Res);
%NN=floor(nx/16);k=0;
NN=floor(nx/16);k=0;
for ii=1:3:nx-NN
    for jj=1:3:ny-NN
        %tmp=gpuArray(abs(recon_Res(jj:jj+NN-1,ii:ii+NN-1,:)));
        tmp=abs(squeeze(recon_Res(jj:jj+NN-1,ii:ii+NN-1,:)));
        bin_tmp = tmp;
        bin_tmp(bin_tmp>0.000000001) = 1;
        s = sum(sum(sum(bin_tmp,1),2),3);
        if s/(NN*NN*nt) > 0.7
            k=k+1;
            Signal(:,k)=squeeze(sum(sum(tmp,1),2));
            %Signal(:,k)= Signal(:,k)/(NN*NN);
            %aux = gpuArray(Signal(:,k));
            temp=abs(fftshift(fft(Signal(:,k))));
            Signal_FFT(:,k) = temp/max(temp(:));
            %Signal_FFT(:,k)=gather(temp);clear temp tmp
        end
    end
end


Res_Peak=squeeze(Signal_FFT(FR_Index,:));
Car_Peak=squeeze(Signal_FFT(FC_Index,:));

ratio_Peak = max(Res_Peak)./max(Car_Peak);
[m,n]=find(Res_Peak==max(Res_Peak(:)));
%n = find(ratio_Peak == max(ratio_Peak));
%m = find(Res_Peak(:,n) == max(Res_Peak(:,n)));

disp(sprintf('Peak requency: %f', Res_Peak(m,n)));

ResFS=F_X(FR_Index);
ResFS=ResFS(m);

disp(sprintf('Respiratory motion frequency: %f', ResFS));
para.ResFS=ResFS;

Res_Signal=Signal(:,n);
Res_Signal_FFT=Signal_FFT(:,n);
Res_Signal=smooth(Res_Signal,7,'moving');

%close all
figure
subplot(2,1,1);plot(time,Res_Signal),title('Respiratory Motion Signal')
subplot(2,1,2);plot(F_X,Res_Signal_FFT),set(gca,'XLim',[-1.5 1.5]),set(gca,'YLim',[-.02 0.08]),
figure,imagescn(abs(recon_Res),[0 .003],[],[],3)


Res_Signal_Long = interp1( linspace(0,1,length(Res_Signal)), Res_Signal, linspace(0,1,para.nt), 'linear');
figure, plot(Res_Signal_Long)
span = double(idivide(int32(para.span),2)*6+1);
Res_Signal_Smooth = smooth(Res_Signal_Long, span, 'lowess');
%figure, plot(Res_Signal_Smooth)
[peak_values,peak_index]= findpeaks(double(Res_Signal_Smooth));
[valley_values,valley_index]= findpeaks(-double(Res_Signal_Smooth));
[peak_values,peak_index,valley_values,valley_index] = SnapExtrema( peak_values,peak_index,valley_values,valley_index, Res_Signal_Long, para.span);

avg_valley = abs(median(valley_values));
max_peak = abs(max(peak_values));

if ResSort  
  Res_Signal1 = InvertRespCurve( Res_Signal_Long, peak_index, valley_index);    
  figure, plot(Res_Signal1)
else
  Res_Signal1 = Res_Signal_Long;
end



% if ResSort
%     for ii=1:length(Res_Signal)-1
%         if Res_Signal(ii)<Res_Signal(ii+1)
%             Res_Signal2(ii)=Res_Signal(ii)*-1;
%         else
%             Res_Signal2(ii)=Res_Signal(ii);
%         end
%     end
%     if Res_Signal(end-1)<Res_Signal(end)
%         Res_Signal2(length(Res_Signal))=Res_Signal(end)*-1;
%     else
%         Res_Signal2(length(Res_Signal))=Res_Signal(end);
%     end
% end


%Res_Signal_new= smooth(Res_Signal,5);

% Res_Signal_new=zeros(nt*4,1);
% Res_Signal_new(1:4:end)=Res_Signal;
% Res_Signal_new(2:4:end)=Res_Signal;
% Res_Signal_new(3:4:end)=Res_Signal;
% Res_Signal_new(4:4:end)=Res_Signal;
