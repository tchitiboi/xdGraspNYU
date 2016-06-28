function  res = MCNUFFT_mskd(k,w,b1,msk)

% Multicoil NUFFT operator
% Based on the NUFFT toolbox from Jeff Fessler
% Input
% k: k-space trajectory
% w: density compensation
% b1: coil sensitivity maps
%
% Ricardo Otazo, 2012
%b1 coil sensitivity map    
% Rebecca Ramb, May 2016 - introduced coil mask

Nd = size(b1(:,:,1));
%trying to weight the different coil elements inside b1
%Nd = size(b1(:,:,[7,8]*0.3));
Jd = [6,6];
Kd = floor([Nd*1.5]);
n_shift = Nd/2;
for tt=1:size(k,3),
	kk=k(:,:,tt);
	om = [real(kk(:)), imag(kk(:))]*2*pi;
	res.st{tt} = nufft_init(om, Nd, Jd, Kd, n_shift,'kaiser');
end
res.adjoint = 0;
res.imSize = size(b1(:,:,1));
res.imSize2 = [size(k,1),size(k,1)];
res.dataSize = size(k);
res.w = sqrt(w);
res.b1 = b1;
res.msk = msk;
res = class(res,'MCNUFFT_mskd');

