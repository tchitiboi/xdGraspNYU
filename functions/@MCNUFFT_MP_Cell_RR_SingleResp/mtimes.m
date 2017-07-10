function ress = mtimes(a,bb)

if a.adjoint,
    % Multicoil non-Cartesian k-space to Cartesian image domain
    % nufft for each coil and time point
    %for resp=1:size(bb,1)        
        for card=1:size(bb,1)
            for label=1:size(bb,2)
                c = bb{card,label};
                aux = a.w{card,label};
                for ch=1:size(c,3)
                    b = c(:,:,ch).*aux;
                    res(:,:,ch) = reshape(nufft_adj(b(:),a.st{card,label})/sqrt(prod(a.imSize2)),a.imSize(1),a.imSize(2));
                end
                ress(:,:,card,label) = sum(res.*conj(a.b1),3)./sum(abs((a.b1)).^2,3);
                clear res
                ress(:,:,card,label) = ress(:,:,card,label).*size(aux,1)*pi/2/size(aux,2);
            end
        end
    %end
else
    % Cartesian image to multicoil non-Cartesian k-space
    ress = cell(size(a.w));
    for card=1:size(bb,3)
       % for resp=1:size(bb,3)
            for label=1:size(bb,4)
                aux = a.w{card,label};
                s = a.dataSize{card,label};
                ress1 = zeros(s(1,1),s(1,2),size(a.b1,3));
                for ch=1:size(a.b1,3)
                    res=bb(:,:,card,label).*a.b1(:,:,ch);
                    ress1(:,:,ch) = reshape(nufft(res,a.st{card,label})/sqrt(prod(a.imSize2)),s(1,1),s(1,2)).*aux;                 
                    clear res
                end
                ress{card,label} = ress1;
                clear ress1
            end
       % end
    end
end