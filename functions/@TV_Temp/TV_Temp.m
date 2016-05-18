function  res = TV_Temp()

%res = TV_Temp()
%
% Implements a temporal finite-differencing operator for dynamic data.
%

res.adjoint = 0;
res = class(res,'TV_Temp');

