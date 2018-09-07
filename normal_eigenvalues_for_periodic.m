function evs = normal_eigenvalues_for_periodic(q,param,tau)
%
% 
    L = param.domain(2) - param.domain(1);
    j = floor(((0:q-1)+1)/2);
    evs = exp(-0.5*(j/tau).^2);
    evs = L/theta3(exp(-(1/(2*tau^2))))*evs;
end