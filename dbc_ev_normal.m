function evs = dbc_ev_normal(param)
%
    a = param.domain(1); b = param.domain(2);
    m = param.m;
    l = param.ev.l;
    v = param.ev.sigma^2;
    
    P = b-a;
    rho = 1/P;        
    
    J = 1:m;
    evs = exp(-2*(pi*J*rho*0.5*l).^2);
    % normalize so that its infinite sum to 1, first.
    evs = evs/(0.5*(-1+theta3(exp(-2*(pi*rho*0.5*l).^2))));

    % obtain max{k(x,x)}=k(c,c) when sum(κ_j)=1.
    c = .5*(a+b);
    Phi = param.ef.fh(c,param);
    kMax = Phi*diag(evs)*Phi';
    
    % k=k(x,x) of the corresponding stationary kernel.
    k = v;
    
    % normalize
    evs = k/kMax*evs;
end