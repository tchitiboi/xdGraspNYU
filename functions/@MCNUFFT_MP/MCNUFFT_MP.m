function  res = MCNUFFT_MP(k,w,b1)

% Multicoil NUFFT operator with multiple respiratory phases
% Based on the NUFFT toolbox from Jeff Fessler
% Input
% k: k-space trajectory
% w: density compensation
% b1: coil sensitivity maps
%
% Ricardo Otazo, 2012
    
Nd = size(b1(:,:,1,1));
Jd = [6,6];
Kd = floor([Nd*1.5]);
n_shift = Nd/2;
for tt=1:size(k,4)
    for resp=1:size(k,3)
        kk=k(:,:,resp,tt);
        om = [real(kk(:)), imag(kk(:))]*2*pi;
        res.st{resp,tt} = nufft_init(om, Nd, Jd, Kd, n_shift,'kaiser');
    end
end
res.adjoint = 0;
res.imSize = size(b1(:,:,1,1));
res.imSize2 = [size(k,1),size(k,1)];
res.dataSize = size(k);
res.w = sqrt(w);
res.b1 = b1;
res = class(res,'MCNUFFT_MP');

