% test wether operator param.E' is the adjoint of operator param.E
% param.E   Operater
% param.y   k-space data
% recon     images space data

p1 = rand(size(recon))+ 1i*rand(size(recon));
p2 = rand(size(param.y)) + 1i*rand(size(param.y));
s1 = param.E*p1;
s2 = param.E'*p2;
p1s2 = p1(:)' * s2(:);
p2s1 = s1(:)' * p2(:);
disp(['p1*s2 = ', num2str(p1s2)]);
disp(['p2*s1 = ', num2str(p2s1)]);
disp(['p1*s2-p2*s1 = ', num2str(p1s2-p2s1)]);