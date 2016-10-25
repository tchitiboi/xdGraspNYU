function ress = mtimes(a,bb)

if a.adjoint,
    % Multicoil non-Cartesian k-space to Cartesian image domain
    % nufft for each coil and time point
    for resp=1:size(bb,1)        
        for card=1:size(bb,2)
            c = bb{resp,card};
            aux = a.w{resp,card};
            for ch=1:size(c,3)
                b = c(:,:,ch).*aux(:,:);
                res(:,:,ch) = reshape(nufft_adj(b(:),a.st{resp,card})/sqrt(prod(a.imSize2)),a.imSize(1),a.imSize(2));
            end
            ress(:,:,resp,card) = sum(res.*conj(a.b1),3)./sum(abs((a.b1)).^2,3);clear res
            ress(:,:,resp,card) = ress(:,:,resp,card).*size(aux,1)*pi/2/size(aux,2);
        end
    end
    %ress=ress.*size(a.w,1)*pi/2/size(a.w,2);
else
    % Cartesian image to multicoil non-Cartesian k-space
    for tt=1:size(bb,4)
        for resp=1:size(bb,3),
            for ch=1:size(a.b1,3),
                res=bb(:,:,resp,tt).*a.b1(:,:,ch);
                ress(:,:,ch,resp,tt) = reshape(nufft(res,a.st{resp,tt})/sqrt(prod(a.imSize2)),a.dataSize(1),a.dataSize(2)).*a.w(:,:,resp,tt);
            end
        end
    end
end