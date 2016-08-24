function [Cardiac_Signal,para]=GetCardiacMotionSignal_HeartBlock(kdata,Traj,DensityComp,b1,nline,para,recon_Car,maskHeart);
%Extract cardiac motion signal from reconstructed low temporal
%resolution images. 
%close all

if ~exist('recon_Car','var')
    nline_car=nline*2;
    nt=floor(size(kdata,2)/nline_car);
    [nx,ny,nc]=size(b1);
    NN=floor(nx/3);%NN=80;
    b1=b1((nx-NN)/2+1:end-(nx-NN)/2,(nx-NN)/2+1:end-(nx-NN)/2,:);
    % if(~bKwic)
    clear kdata_Under Traj_Under DensityComp_Under
    for ii=1:nt
        kdata_Under(:,:,:,ii)=kdata(:,(ii-1)*nline_car+1:ii*nline_car,:);
        Traj_Under(:,:,ii)=Traj(:,(ii-1)*nline_car+1:ii*nline_car);
        DensityComp_Under(:,:,ii)=DensityComp(:,(ii-1)*nline_car+1:ii*nline_car);
    end
    E=MCNUFFT(Traj_Under,DensityComp_Under,b1);
    [nx,ntviews,nc,nt]=size(kdata_Under);
    recon_Car=E'*double(kdata_Under.*repmat(kaiser(nx,20),[1,ntviews,nc,nt]));
    clear Traj_Under DensityComp_Under E
    % else
    %     [recon_Car,kdata_Under,Traj_Under,DensityComp_Under,kwicmask,kwicdcf] = apply_kwic(kdata,Traj,DensityComp,b1,nline_car,0);
    % end
end

[nx,ny,nt]=size(recon_Car);

para.LF_H=0.8;para.HF_H=2;%%% initial heart rate range
[HF_Index, F_X] = selectCardiacMotionFrequencies(para, nt);

time_series = abs(recon_Car);
for t = 1:nt
 %sigmoid filtering
 tmp(:,:,t)= repmat(1,[nx,ny])./(1 + exp(-(time_series(:,:,t)-repmat(0.8,[nx, ny]))/0.3));
end

if ~exist('maskHeart','var')
maskHeart = tmc_localizeHeart(time_series, HF_Index);
end
para.maskHeart = maskHeart;

tmp = time_series.*repmat(maskHeart,[1 1 nt]);
tmp = tmp/max(max(max(tmp)));

figure,imagescn(abs(tmp),[0 .6*max(tmp(:))],[],[],3)

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
figure,imagescn(abs(recon_Car),[0 .6*max(abs(recon_Car(:)))],[],[],3)