function evs = studentT_eigenvalues_for_dirichlet(param)
% Since simple expression of sum of {t_\nu(j|\sigma} is not known,
% normalization is performed based on sum(evs).

    a = param.domain(1); b = param.domain(2);
    m = param.m;
    
    v = param.ev.sigma^2;
    scale = param.ev.scale;
    nu = param.ev.nu;
    
    J = 1:m;

    % normalize so that sum to 1, first.
    t = makedist('tLocationScale', ...
        'mu',0,'sigma',scale,'nu',nu);
    evs = pdf(t,J);
    evs = evs/sum(evs);
    
    % obatain max{k(x,x)}=k(c,c) when sum(κ_j)=1.
    c = .5*(a+b);
    Phi = param.ef.fh(c,param);
    kMax = Phi*diag(evs)*Phi';
    
    % k(x,x) of the corresponding stationary kernel.
    k = v;
    
    % normalize
    evs = k/kMax*evs;
end