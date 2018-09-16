function evs = pbc_ev_normal(param)
% pbc_ev_normal returns eigenvalues defined using Normal distribution
% for the eigensystem under the periodic boundary condition.
% 
    a = param.domain(1); b = param.domain(2);
    m = param.m;
    v = param.ev.sigma^2;
    l = param.ev.l;
    P = b-a;
    rho = 1/P;
    
    S = floor(((0:m-1)+1)/2);
    evs = v*exp(-2*(pi*S*rho*l).^2);
    evs = P/theta3(exp(-2*(pi*rho*l).^2))*evs;
end