function [Cardiac_Signal,para]=GetCardiacMotionSignal_HeartBlock(kdata,Traj,DensityComp,b1,nline,para);
%Extract cardiac motion signal from reconstructed low temporal
%resolution images. 
close all
nline_car=nline*2;
nt=floor(size(kdata,2)/nline_car);
for ii=1:nt
    kdata_Under(:,:,:,ii)=kdata(:,(ii-1)*nline_car+1:ii*nline_car,:);
    Traj_Under(:,:,ii)=Traj(:,(ii-1)*nline_car+1:ii*nline_car);
    DensityComp_Under(:,:,ii)=DensityComp(:,(ii-1)*nline_car+1:ii*nline_car);
end
[nx,ny,nc]=size(b1);
NN=floor(nx/3);%NN=80;
b1=b1((nx-NN)/2+1:end-(nx-NN)/2,(nx-NN)/2+1:end-(nx-NN)/2,:); 
%b1 sensitivity map for all coils
%code to try to give weight in the reconstruction to the different coils
%b1=b1((nx-NN)/2+1:end-(nx-NN)/2,(nx-NN)/2+1:end-(nx-NN)/2,[1,3])*1;
%b1=b1((nx-NN)/2+1:end-(nx-NN)/2,(nx-NN)/2+1:end-(nx-NN)/2,[4,5,6,9])*0.5;
%b1=b1((nx-NN)/2+1:end-(nx-NN)/2,(nx-NN)/2+1:end-(nx-NN)/2,[2,7,8]*0.01,[1,3,4,5,6,9,10])*0.1;
%b1=b1((nx-NN)/2+1:end-(nx-NN)/2,(nx-NN)/2+1:end-(nx-NN)/2,[2,7,8]*1);
%b1=b1((nx-NN)/2+1:end-(nx-NN)/2,(nx-NN)/2+1:end-(nx-NN)/2,[1,3,4,5,6,9,10])*0.01;
E=MCNUFFT(Traj_Under,DensityComp_Under,b1);
[nx,ntviews,nc,nt]=size(kdata_Under);
recon_Car=E'*double(kdata_Under.*repmat(kaiser(nx,20),[1,ntviews,nc,nt]));
clear Traj_Under DensityComp_Under E

[nx,ny,nt]=size(recon_Car);

[HF_Index, F_X] = selectCardiacMotionFrequencies(para, nt);

time_series = abs(recon_Car);
maskHeart = tmc_localizeHeart(time_series, HF_Index);

tmp = time_series.*repmat(maskHeart,[1 1 nt]);

%tmp = medfilt1(tmp, 3, [], 3);

tmp = tmp/max(max(max(tmp)));

for t = 1:nt
 %sigmoid filtering
 tmp(:,:,t)= repmat(1,[nx,ny])./(1 + exp(-(tmp(:,:,t)-repmat(0.65,[nx, ny]))/0.2));
end

figure,imagescn(abs(tmp),[0 1],[],[],3)

Signal=squeeze(sum(sum(tmp,1),2));
temp=abs(fftshift(fft(Signal)));
Signal_FFT=temp/max(temp(:));clear temp tmp
% NN=70;k=0;
% for ii=1:1:nx-NN
%     for jj=1:1:ny-NN
%         k=k+1;
%         tmp=abs(recon_Res(jj:jj+NN-1,ii:ii+NN-1,:));
%         Signal(:,k)=squeeze(sum(sum(tmp,1),2));
%         temp=abs(fftshift(fft(Signal(:,k))));
%         Signal_FFT(:,k)=temp/max(temp(:));clear temp tmp
%     end
% end


Heart_Peak=squeeze(Signal_FFT(HF_Index,:));

[m,n]=find(Heart_Peak==max(Heart_Peak(:)));
HeartFS=F_X(HF_Index);
HeartFS=HeartFS(m);
% 
disp(sprintf('Cardiac motion frequency: %f', HeartFS));
para.HeartFS=HeartFS;

Cardiac_Signal=Signal;
Cardiac_Signal_FFT=Signal_FFT;
Cardiac_Signal=smooth(Cardiac_Signal,4,'lowess');

Fs=1/para.TR;
fc_high = para.ResFS*1.5;
Filter_Bandpass = fir1(size(Cardiac_Signal,1), 2*fc_high/Fs,'low');
Signal_Filtered=conv([Cardiac_Signal;Cardiac_Signal;Cardiac_Signal],Filter_Bandpass);
temp = Signal_Filtered(size(Cardiac_Signal,1)*1.5+1:size(Cardiac_Signal,1)*2.5,:);
Cardiac_Signal=Cardiac_Signal/max(Cardiac_Signal);
temp=temp./max(temp);
Cardiac_Signal=Cardiac_Signal./temp;clear temp
Cardiac_Signal=imresize(Cardiac_Signal,[nt*2,1]);

time = para.TR:para.TR:nt*2*para.TR;

close all
figure
subplot(2,1,1);plot(time,Cardiac_Signal),title('Cardiac Motion Signal')
subplot(2,1,2);plot(F_X,Cardiac_Signal_FFT),set(gca,'XLim',[-2 2]),set(gca,'YLim',[-.02 0.08]),
figure,imagescn(abs(recon_Car),[0 .001],[],[],3)
