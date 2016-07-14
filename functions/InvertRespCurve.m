function [ Res_Signal1 ] = InvertRespCurve( Res_Signal, peak_index, valley_index)

vIndex = 1;

if valley_index(vIndex) < peak_index(vIndex)
    for i = 1 : valley_index(vIndex)
        if Res_Signal(i)<Res_Signal(i+1)
            Res_Signal1(i)=Res_Signal(i)*-1;
        else
            Res_Signal1(i)=Res_Signal(i);
        end
    end
    for i = valley_index(vIndex)+1 : peak_index(vIndex)-1
        Res_Signal1(i) = Res_Signal(i)*-1;
    end
else
    for i = 1 : peak_index(vIndex)-1
        Res_Signal1(i) = Res_Signal(i)*-1;
    end
end

for ii = 1 : length(peak_index) - 1
    while valley_index(vIndex) < peak_index(ii)
        vIndex = vIndex + 1;
    end
    for i = peak_index(ii) : valley_index(vIndex)
        Res_Signal1(i) = Res_Signal(i);
    end
    for i = valley_index(vIndex)+1 : peak_index(ii+1)-1
        Res_Signal1(i) = Res_Signal(i)*-1;
    end
end

ii = length(peak_index);
vIndex = length(valley_index);
if valley_index(vIndex) > peak_index(ii)
    for i = peak_index(ii) : valley_index(vIndex)
        Res_Signal1(i) = Res_Signal(i);
    end
    for i = valley_index(vIndex)+1 : length(Res_Signal)
        Res_Signal1(i) = Res_Signal(i)*-1;
    end
else
    for i = peak_index(ii) : length(Res_Signal)
        Res_Signal1(i) = Res_Signal(i);
    end
end

Res_Signal1 = Res_Signal1 - min(Res_Signal1);
Res_Signal1 = permute(Res_Signal1, [2,1]);

end

