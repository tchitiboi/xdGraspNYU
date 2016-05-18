function para=TimeStamp(para);
%Generate the time and frequency stamp

para.time = para.TR:para.TR:para.nt*para.TR;
F_S = 1/para.TR;F_X = 0:F_S/(para.nt-1):F_S;
F_X=F_X-F_S/2;  %%% frequency after FFT of the motion signal
if mod(para.nt,2)==0
    F_X=F_X+F_X(para.nt/2);
end
para.F_X=F_X;
