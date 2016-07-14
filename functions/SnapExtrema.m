function [peak_values,peak_index,valley_values,valley_index] = SnapExtrema( peak_values,peak_index,valley_values,valley_index, Res_Signal, span)
% Shifts the extrema determined on the smoothed signal to the positions on
% the original signal.

for i = 1 : length(valley_index)
   for ii = max(valley_index(i) - span, 1) : min(valley_index(i) + span, length(Res_Signal))
      if (Res_Signal(ii) < valley_values(i))
         valley_values(i) =  Res_Signal(ii);
         valley_index(i) == ii;
      end
   end
end

for i = 1 : length(peak_index)
   for ii = max(peak_index(i) - span, 1) : min(peak_index(i) + span, length(Res_Signal))
      if Res_Signal(ii) > peak_values(i)
         peak_values(i) =  Res_Signal(ii);
         peak_index(i) == ii;
      end
   end
end

end

