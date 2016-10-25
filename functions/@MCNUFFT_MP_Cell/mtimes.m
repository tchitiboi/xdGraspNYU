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
    ress = cell(size(a.w));
    for card=1:size(bb,4)
        for resp=1:size(bb,3)
            aux = a.w{resp,card};
            s = a.dataSize{resp,card};
            %c = bb{resp,card};
            for ch=1:size(a.b1,3)
                res=bb(:,:,resp,card).*a.b1(:,:,ch);
                ress1(:,:,ch) = reshape(nufft(res,a.st{resp,card})/sqrt(prod(a.imSize2)),s(1,1),s(1,2)).*aux;
                ress{resp,card} = ress1;
                clear res ress1
            end
        end
    end
end