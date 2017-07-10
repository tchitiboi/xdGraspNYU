function [x] = CSL1NlCg_Cell_w_GPU(x0,param)
% 
% res = CSL1NlCg(param)

% starting point
x=x0;

% line search parameters
maxlsiter = 6;
gradToll = 1e-8 ;
param.l1Smooth = 1e-15;	
alpha = 0.01;  
beta = 0.6;
t0 = 1 ; 
k = 0;
% compute g0  = grad(f(x))
g0 = grad(x,param);

%  %% Preparation
%     fprintf(' Performing gradient check \n')
%     
%     %x = x0;
%     rng(1);
%     dx = rand(size(x));
%     % dx = rand(size(x)) + 1i*rand(size(x));
%     % h = rand(size(x)) + 1i*rand(size(x));
%     t = 1e-10;
%     %t=t0;
%     
%     % Calculate variables
%     % res_F       = objective_fw_lin_grasp_fm(x, dx, 0, kdata, par);
%     %res_F_plus  = objective_fw_lin_grasp_fm(x + t*dx, kdata, par);
%     %res_F_minus = objective_fw_lin_grasp_fm(x - t*dx, kdata, par);
%     
% %     res_F       = objective(x, dx, 0, param);
%     res_F_plus  = objective(x, dx, t, param);
%     res_F_minus = objective(x, dx, -t, param);
%     
% %     J_fd        = (res_F_plus - res_F)/(t); % forward
%     % J_fd        = (res_F - res_F_minus)/(t); % backward
%     J_fd        = (res_F_plus - res_F_minus)/(2*t); % central
%     
%     g = grad(x, param);
%     J_an = g(:)'*dx(:);
%     
%     J_an = real(J_an);
%     ratio = J_fd/J_an;
%     
%     disp([' Finite difference = ', num2str(J_fd), ', Analytical = ', num2str(J_an), ', ratio = ', num2str(ratio)]);
%     fprintf(' ------------------------- \n')

dx = -g0;
% iterations
while(1)

    % backtracking line-search
	f0 = objective(x,dx,0,param);
	t = t0;
    f1 = objective(x,dx,t,param);
	lsiter = 0;
	while (f1 > f0 - alpha*t*abs(g0(:)'*dx(:)))^2 & (lsiter<maxlsiter)
		lsiter = lsiter + 1;
		t = t * beta;
		f1 = objective(x,dx,t,param);
	end

	% control the number of line searches by adapting the initial step search
	if lsiter > 2, t0 = t0 * beta;end 
	if lsiter<1, t0 = t0 / beta; end

	x = (x + t*dx);

    % print some numbers for debug purposes	
    disp(sprintf('%d   , obj: %f, L-S: %d', k,f1,lsiter));

    %conjugate gradient calculation
	g1 = grad(x,param);
	bk = g1(:)'*g1(:)/(g0(:)'*g0(:)+eps);
	g0 = g1;
	dx =  - g1 + bk* dx;
	k = k + 1;

	% stopping criteria (to be improved)
	if (k > param.nite) || (norm(dx(:)) < gradToll), break;end

end
return;

function res = objective(x,dx,t,param) %**********************************

% L2-norm part
L2Obj=0;
aux = x+t*dx;
for resp = 1:size(param.E1,1)
    for card = 1:size(param.E1,2)
       p(1,:,1) = param.SGW{resp,card};       
       aux1 = (param.E1{resp,card}*squeeze(aux(:,:,resp,card)));%*0.2081;%.*sqrt(param.nx*pi/2/size(param.y{resp,card},2))/param.b1_scalar;
       w = bsxfun(@times,p, reshape(aux1, [param.nx, size(aux1,1)/param.nx, size(aux1,2)]) - param.y{resp,card});
       L2Obj = L2Obj + w(:)'*w(:);
       clear p aux1 w
    end
end

if param.TVWeight
    w = param.TV*(x+t*dx); 
    TVObj = sum((w(:).*conj(w(:))+param.l1Smooth).^(1/2));
else
    TVObj = 0;
end

if param.L1Weight
    w = param.W*(x+t*dx); 
    WObj = sum((w(:).*conj(w(:))+param.l1Smooth).^(1/2));
else
    WObj = 0;
end

res=L2Obj+param.TVWeight*TVObj+param.L1Weight*WObj;

function g = grad(x,param)%***********************************************

% L2-norm part
L2Grad = zeros(size(x));

for resp = 1:size(param.E1,1)
    for card = 1:size(param.E1,2)
        p(1,:,1) = param.SGW{resp,card}.^2;
        aux1 = param.E1{resp,card}*x(:,:,resp,card);%*0.2081;%.*sqrt(param.nx*pi/2/size(param.y{resp,card},2))/param.b1_scalar;
        aux = bsxfun(@times,reshape(aux1, [param.nx, size(aux1,1)/param.nx, size(aux1,2)]) - param.y{resp,card},p);
        L2Grad(:,:,resp,card) = 2.*(param.E1{resp,card}'*double(reshape(aux, [size(aux,1)*size(aux,2),size(aux,3)])));
        L2Grad(:,:,resp,card) = L2Grad(:,:,resp,card)./sum(abs((param.b1)).^2,3).*param.nx*pi/2/size(aux,2);%*0.2081;
        clear p aux aux1
    end
end

% L2-norm part
% L2Grad = 2.*(param.E'*(param.E*x-param.y));

if param.TVWeight
    w = param.TV*x;
    TVGrad = param.TV'*(w.*(w.*conj(w)+param.l1Smooth).^(-0.5));
else
    TVGrad=0;
end

if param.L1Weight
    w = param.W*x;
    WGrad = param.W'*(w.*(w.*conj(w)+param.l1Smooth).^(-0.5));
else
    WGrad=0;
end

g=L2Grad+param.TVWeight*TVGrad+param.L1Weight*WGrad;
