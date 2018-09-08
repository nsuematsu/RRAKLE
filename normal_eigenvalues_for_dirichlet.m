function evs = normal_eigenvalues_for_dirichlet(q,param,tau)
%
% 
    L = param.domain(2) - param.domain(1);
    J = 1:q;
    evs = exp(-0.5*(J/tau).^2);
    evs = L/(0.5*(theta3(exp(-(1/(2*tau^2))))-1))*evs;
end