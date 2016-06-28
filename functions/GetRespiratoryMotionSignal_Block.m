function [Res_Signal_new,para]=GetRespiratoryMotionSignal_Block(kdata,Traj,DensityComp,b1,nline,para,ResSort);

global osf % oversampling: 1.5 1.25
global wg % kernel width: 5 7
global sw % parallel sectors' width: 12 16

%Extract respiratory motion signal from reconstructed low temporal
%resolution images. 
%close all
nline_res=nline*4;
nt=floor(size(kdata,2)/nline_res);
[nx,ny,nc]=size(b1);
imwidth = nx;

for ii=1:nt
    kdata_Under(:,:,:,ii)=kdata(:,(ii-1)*nline_res+1:ii*nline_res,:);
    Traj_Under(:,:,ii)=Traj(:,(ii-1)*nline_res+1:ii*nline_res);
    DensityComp_Under(:,:,ii)=DensityComp(:,(ii-1)*nline_res+1:ii*nline_res);
end

[nx,ntviews,nc,nt]=size(kdata_Under);
kdata_Under = kdata_Under.*repmat(kaiser(nx,20),[1,ntviews,nc,nt]);

if para.flag 
    tic
    NN=192; %?????
    b1=b1((nx-NN)/2+1:end-(nx-NN)/2,(nx-NN)/2+1:end-(nx-NN)/2,:);

    E1=MCNUFFT(Traj_Under,DensityComp_Under,b1);
    [nx,ntviews,nc,nt]=size(kdata_Under);
    recon_Res=E1'*double(kdata_Under);
    timeMin = toc;
    timeMin = timeMin/60
    clear Traj_Under DensityComp_Under E
else    
    tic
    for ii=1:nt
        E = gpuNUFFT([real(col(Traj_Under(:,:,ii))), imag(col(Traj_Under(:,:,ii)))]',col(DensityComp_Under(:,:,ii)),osf,wg,sw,[imwidth,imwidth],b1,true);

        kdata_Under_Res(:,:) = reshape(kdata_Under(:,:,:,ii),[size(kdata_Under,1)*size(kdata_Under,2),nc]);
        rec_Res = E'*double(kdata_Under_Res);

        recon_Res(:,:,ii) = rec_Res;   
    end

    time = toc;
    time = time/60
end

figure,imagescn(abs(recon_Res),[0 .00003],[],[],3)


TR=para.TR*4;
time = TR:TR:nt*TR;
F_S = 1/TR;F_X = 0:F_S/(nt-1):F_S;
F_X=F_X-F_S/2;  %%% frequency after FFT of the motion signal
if mod(nt,2)==0
    F_X=F_X+F_X(nt/2);
end
FR_Index=find(F_X<para.HF_R & F_X>para.LF_R);
FC_Index=find(F_X<para.HF_H & F_X>para.LF_H);

if(~para.flag)
    g = gpuDevice(1);
    reset(g);
end

[nx,ny,nt]=size(recon_Res);
NN=floor(nx/16);k=0;
for ii=1:3:nx-NN
    for jj=1:3:ny-NN
        k=k+1;
        %tmp=gpuArray(abs(recon_Res(jj:jj+NN-1,ii:ii+NN-1,:)));
        tmp=abs(recon_Res(jj:jj+NN-1,ii:ii+NN-1,:));
        Signal(:,k)=squeeze(sum(sum(tmp,1),2));
        %aux = gpuArray(Signal(:,k));
        temp=abs(fftshift(fft(Signal(:,k))));
        Signal_FFT(:,k) = temp/max(temp(:));
        %Signal_FFT(:,k)=gather(temp);clear temp tmp
    end
end


Res_Peak=squeeze(Signal_FFT(FR_Index,:));
Car_Peak=squeeze(Signal_FFT(FC_Index,:));

ratio_Peak = max(Res_Peak)./max(Car_Peak);
n = find(ratio_Peak == max(ratio_Peak));
m = find(Res_Peak(:,n) == max(Res_Peak(:,n)));

%[m,n]=find(Res_Peak==max(Res_Peak(:)));
ResFS=F_X(FR_Index);
ResFS=ResFS(m);


disp(sprintf('Respiratory motion frequency: %f', ResFS));
para.ResFS=ResFS;

Res_Signal=Signal(:,n);
Res_Signal_FFT=Signal_FFT(:,n);
Res_Signal=smooth(Res_Signal,5,'moving');

%%%%%%%%%%%%%%%%%%%%%%%NOT USE
% % % %Band-pass filtering to the respiratory motion signal
% % % Fs=1/TR;
% % % fc_low = ResFS/3.5;
% % % fc_high = ResFS*2.5;
% % % Filter_Bandpass = fir1(size(Res_Signal,1), [2*fc_low/Fs 2*fc_high/Fs]);
% % % Signal_Filtered=conv([Res_Signal;Res_Signal;Res_Signal],Filter_Bandpass);
% % % Res_Signal_new = Signal_Filtered(size(Res_Signal,1)*1.5+1:size(Res_Signal,1)*2.5,:);
%%%%%%%%%%%%%%%%%%%%%%%%%

close all
figure
subplot(2,1,1);plot(time,Res_Signal),title('Respiratory Motion Signal')
subplot(2,1,2);plot(F_X,Res_Signal_FFT),set(gca,'XLim',[-1.5 1.5]),set(gca,'YLim',[-.02 0.08]),
figure,imagescn(abs(recon_Res),[0 .003],[],[],3)

if ResSort
    for ii=1:length(Res_Signal)-1
        if Res_Signal(ii)<Res_Signal(ii+1)
            Res_Signal1(ii)=Res_Signal(ii)*-1;
        else
            Res_Signal1(ii)=Res_Signal(ii);
        end
    end
    if Res_Signal(end-1)<Res_Signal(end)
        Res_Signal1(length(Res_Signal))=Res_Signal(end)*-1;
    else
        Res_Signal1(length(Res_Signal))=Res_Signal(end);
    end
    Res_Signal=Res_Signal1;clear Res_Signal1
end
Res_Signal_new=imresize(Res_Signal,[nt*4,1]);
% Res_Signal_new=zeros(nt*4,1);
% Res_Signal_new(1:4:end)=Res_Signal;
% Res_Signal_new(2:4:end)=Res_Signal;
% Res_Signal_new(3:4:end)=Res_Signal;
% Res_Signal_new(4:4:end)=Res_Signal;
