function [nrm] = power_it_cell(K,Kt,dims,n)

x = rand(dims);

for i=1:n
    aux = K(x);
    x = Kt(aux);
    nrm = sqrt(sum(abs(x(:)).^2));
    x = x./nrm;
    if rem(i,10)==0
        display(['At it: ',num2str(i)]);
    end
end

%Nrm converges to |K|^2, hence:

nrm = sqrt(nrm);

