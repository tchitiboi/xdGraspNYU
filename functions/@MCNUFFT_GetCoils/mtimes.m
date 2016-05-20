function ress = mtimes(a,bb)

if a.adjoint,
    % Multicoil non-Cartesian k-space to Cartesian image domain
    % nufft for each coil and time point
    for tt=1:size(bb,4),
        for ch=1:size(bb,3),
            b = bb(:,:,ch,tt).*a.w(:,:,tt);
            %res(:,:,ch) = reshape(nufft_adj(b(:),a.st{tt})/sqrt(prod(a.imSize2)),a.imSize(1),a.imSize(2));
            res(:,:,ch) = nufft_adj(b(:),a.st{tt})/sqrt(prod(a.imSize2));
            res(:,:,ch) = reshape(res(:,:,ch), a.imSize(1), a.imSize(2));%*conj(a.b1(:,:,ch))./(abs(a.b1(:,:,ch)).^2);
        end
        %ress(:,:,tt)=sum(res.*conj(a.b1),3)./sum(abs((squeeze(a.b1))).^2,3);
        ress(:,:,:,tt)=res.*conj(a.b1);%./(abs(a.b1).^2);
        %clear res
    end
    ress=ress.*size(a.w,1)*pi/2/size(a.w,2);
else
    % Cartesian image to multicoil non-Cartesian k-space
    for tt=1:size(bb,4),
        for ch=1:size(a.b1,3),
            res = bb(:,:,ch,tt).*a.b1(:,:,ch); 
            res = nufft(res,a.st{tt})/sqrt(prod(a.imSize2));
            res = reshape(res,a.dataSize(1),a.dataSize(2));
            ress(:,:,ch,tt) = res.*a.w(:,:,tt);
        end
    end
end

