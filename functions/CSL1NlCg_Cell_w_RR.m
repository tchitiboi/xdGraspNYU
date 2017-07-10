function [x] = CSL1NlCg_Cell_w_RR(x0,param,weightMat)
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
g0 = grad(x,param,weightMat);
dx = -g0;

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

% iterations
while(1)

    % backtracking line-search
	f0 = objective(x,dx,0,param,weightMat);
	t = t0;
    f1 = objective(x,dx,t,param,weightMat);
	lsiter = 0;
	while (f1 > f0 - alpha*t*abs(g0(:)'*dx(:)))^2 & (lsiter<maxlsiter)
		lsiter = lsiter + 1;
		t = t * beta;
		f1 = objective(x,dx,t,param,weightMat);
	end

	% control the number of line searches by adapting the initial step search
	if lsiter > 2, t0 = t0 * beta;end 
	if lsiter<1, t0 = t0 / beta; end

	x = (x + t*dx);

    % print some numbers for debug purposes	
    disp(sprintf('%d   , obj: %f, L-S: %d', k,f1,lsiter));

    %conjugate gradient calculation
	g1 = grad(x,param,weightMat);
	bk = g1(:)'*g1(:)/(g0(:)'*g0(:)+eps);
	g0 = g1;
	dx =  - g1 + bk* dx;
	k = k + 1;

	% stopping criteria (to be improved)
	if (k > param.nite) || (norm(dx(:)) < gradToll), break;end

end
return;

function res = objective(x,dx,t,param,weightMat) %**********************************

% L2-norm part
aux = param.E*(x+t*dx);
w = cell(size(aux));
L2Obj=0;
for resp = 1:size(aux,1)
    for card = 1:size(aux,2)
        for label = 1:size(aux,3)
           p(1,:,1) = param.SGW{resp,card,label};
           w{resp,card,label} = bsxfun(@times,p,aux{resp,card,label} - param.y{resp,card,label});
           L2Obj = L2Obj + w{resp,card,label}(:)'*w{resp,card,label}(:);
           clear p
        end
    end
end
%L2Obj=w(:)'*w(:);
clear aux w

TVObj = 0;
if param.TVWeight
    w = param.TV*(x+t*dx); 
    for resp = 1:size(w,3)
        for card = 1:size(w,4)
            for label = 1:size(w,5)
                TVObj = sum(sum(sum((w(:,:,resp,card,label).*conj(w(:,:,resp,card,label))...
                    +param.l1Smooth).^(1/2)*weightMat(resp,card,label))));
            end
        end
    end
end

WObj = 0;
if param.L1Weight
    w = param.W*(x+t*dx); 
    for resp = 1:size(w,3)
        for card = 1:size(w,4)
            for label = 1:size(w,5)
                WObj = sum(sum(sum((w(:,:,resp,card,label).*conj(w(:,:,resp,card,label))...
                    +param.l1Smooth).^(1/2)*weightMat(resp,card,label))));
            end
        end
    end
end

W1Obj = 0;
if param.L1Weight1
    w = param.W1*(x+t*dx); 
    for resp = 1:size(w,3)
        for card = 1:size(w,4)
            for label = 1:size(w,5)
                W1Obj = sum(sum(sum((w(:,:,resp,card,label).*conj(w(:,:,resp,card,label))...
                    +param.l1Smooth).^(1/2)*weightMat(resp,card,label))));
            end
        end
    end
end

res=L2Obj+param.TVWeight*TVObj+param.L1Weight*WObj+param.L1Weight1*W1Obj;

function g = grad(x,param,weightMat)%***********************************************

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

if param.L1Weight1
    w = param.W1*x;
    W1Grad = param.W1'*(w.*(w.*conj(w)+param.l1Smooth).^(-0.5));
else
    W1Grad=0;
end

% extended weights
weightMatExt = zeros(size(WGrad));
% L2-norm part
aux = param.E*x;
aux1 = cell(size(aux));
for resp = 1:size(aux,1)
    for card = 1:size(aux,2)
        for label = 1:size(aux,3)
            p(1,:,1) = param.SGW{resp,card,label}.^2;
            aux1{resp,card,label} = bsxfun(@times,aux{resp,card,label}-param.y{resp,card,label},p);
            weightMatExt(:,:,resp,card,label) = weightMat(resp,card,label);
            clear p
        end
    end
end
L2Grad = 2.*(param.E'*aux1);
clear aux aux1

g=L2Grad+(param.TVWeight*TVGrad+param.L1Weight*WGrad+param.L1Weight1*W1Grad).*weightMatExt;
