function  res = TV_Temp3D()

%res = TV_Temp3D()
%
% Implements a spatial finite-differencing operator for dynamic data.
%
% Ricardo Otazo 2008

res.adjoint = 0;
res = class(res,'TV_Temp3D');

