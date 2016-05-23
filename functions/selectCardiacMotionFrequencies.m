function [FH_Index, F_X] = selectCardiacMotionFrequencies(para, nt)

TR=para.TR*2;
time = TR:TR:nt*TR;
F_S = 1/TR;F_X = 0:F_S/(nt-1):F_S;
F_X=F_X-F_S/2;  %%% frequency after FFT of the motion signal
if mod(nt,2)==0
    F_X=F_X+F_X(nt/2);
end
% 
FH_Index=find(F_X<para.HF_H & F_X>para.LF_H);

end

