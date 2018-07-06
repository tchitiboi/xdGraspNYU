function [Cardiac_Signal,para,maskHeart]=GetCardiacMotionSignal_HeartBlock(kdata,Traj,DensityComp,b1,nline,para,recon_Car);
%Extract cardiac motion signal from reconstructed low temporal
%resolution images. 
close all

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

para.LF_H=1;para.HF_H=1.7;%%% initial heart rate range
[HF_Index, F_X] = selectCardiacMotionFrequencies(para, nt);

time_series = abs(recon_Car);
for t = 1:nt
 %sigmoid filtering
 tmp(:,:,t)= repmat(1,[nx,ny])./(1 + exp(-(time_series(:,:,t)-repmat(0.8,[nx, ny]))/0.3));
end

% tmp1 = zeros(size(tmp));
% for t = 2:(size(tmp,3)-1)
%   tmp1(:,:,t) = (tmp(:,:,t-1)+3*tmp(:,:,t)+tmp(:,:,t+1))/5;
% end
% 
% tmp1(:,:,1) = (3*tmp(:,:,1)+tmp(:,:,2))/4;
% tmp1(:,:,size(tmp,3)) = (tmp(:,:,size(tmp,3)-1)+tmp(:,:,size(tmp,3)))/4;
% tmp = tmp1;

% tmp = medfilt1(tmp,3,[],3);
% time_series = tmp;

maskHeart = tmc_localizeHeart(tmp, HF_Index);
k=0;
for border_size = 0:1:50
    k=k+1;
    border_img = maskHeart;
        for x = 1:size(maskHeart,1)
          for y = 1:size(maskHeart,2)
              if (x < border_size || x > size(maskHeart,1)-border_size) ...
                 || (y < border_size || y > size(maskHeart,1)-border_size)
                  border_img(x,y) = 0;
              else
                  border_img(x,y) = 1;  
              end
          end
        end
    new_img = maskHeart.*border_img;
    tmp1 = time_series.*repmat(new_img,[1 1 nt]);
    tmp1 = tmp1/max(max(max(tmp1)));
    
    %indeces = find(tmp1>0.85);
    %tmp1(indeces) = 0;

    Signal(:,k)=squeeze(sum(sum(tmp1,1),2));
    temp=abs(fftshift(fft(Signal(:,k))));
    Signal_FFT(:,k)=temp/max(temp(:));

end

figure,imagescn(abs(time_series),[0 .6*max(tmp(:))],[],[],3)


% time_series = tmp;
% clear temp tmp
% 
% NN=30;k=0;
% for ii=1:5:nx-NN
%     for jj=1:5:ny-NN
%         k=k+1;
%         tmp=abs(time_series(jj:jj+NN-1,ii:ii+NN-1,:));
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
disp(sprintf('Peak: %f', max(Heart_Peak(:))));
para.HeartFS=HeartFS;

% Cardiac_Signal=Signal(:,n);
% Cardiac_Signal_FFT=Signal_FFT(:,n);

Cardiac_Signal=Signal(:,n);
Cardiac_Signal_FFT=Signal_FFT(:,n);
Cardiac_Signal=smooth(Cardiac_Signal,4,'lowess');

Fs=1/para.TR;
% fc_high = para.ResFS*1.5;
% Filter_Bandpass = fir1(size(Cardiac_Signal,1), 2*fc_high/Fs,'low');
% Signal_Filtered=conv([Cardiac_Signal;Cardiac_Signal;Cardiac_Signal],Filter_Bandpass);
% temp = Signal_Filtered(size(Cardiac_Signal,1)*1.5+1:size(Cardiac_Signal,1)*2.5,:);
Cardiac_Signal=Cardiac_Signal/max(Cardiac_Signal);
% temp=temp./max(temp);
% Cardiac_Signal=Cardiac_Signal./temp;clear temp
Cardiac_Signal=imresize(Cardiac_Signal,[nt*2,1]);

time = para.TR:para.TR:nt*2*para.TR;

%close all
figure
subplot(2,1,1);plot(time,Cardiac_Signal),title('Cardiac Motion Signal')
subplot(2,1,2);plot(F_X,Cardiac_Signal_FFT),set(gca,'XLim',[-2 2]),set(gca,'YLim',[-.02 0.08]),
figure,imagescn(abs(recon_Car),[0 .6*max(abs(recon_Car(:)))],[],[],3)