function  res = MCNUFFT_MP_Cell_RR(k, w, b1)

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
for resp=1:size(k,1)    
    for card=1:size(k,2)    
        for label=1:size(k,3)    
            kk = k{resp,card,label};       
            om = [real(kk(:)), imag(kk(:))]*2*pi;
            res.st{resp,card,label} = nufft_init(om, Nd, Jd, Kd, n_shift,'kaiser');
        end
    end
end
res.adjoint = 0;
res.imSize = size(b1(:,:,1,1));
res.imSize2 = [size(k{1,1},1),size(k{1,1},1)];
aux = cell(size(w));
res.dataSize = cell(size(k));
for resp=1:size(w,1)
    for card=1:size(w,2)
        for label=1:size(w,3)
            mat = w{resp,card,label};
            aux{resp,card,label} = sqrt(mat);
            s = size(k{resp,card,label});
            res.dataSize{resp,card,label} = s;
        end
    end
end
res.w = aux;
res.b1 = b1;
res = class(res,'MCNUFFT_MP_Cell_RR');

