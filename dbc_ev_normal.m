function [evs,gevs] = dbc_ev_normal(param)
%
    a = param.domain(1); b = param.domain(2);
    m = param.m;
    l = param.ev.l;
    v = param.ev.sigma^2;
    
    P = b-a;
    rho = 1/P;        
    
    F = .5*(1:m);
    e = v*exp(-2*(pi*F*rho*l).^2);

    % obtain max{k(x,x)}=k(c,c)
    c = .5*(a+b);
    Phi = param.ef.fh(c,param);
    kMax = Phi*diag(e)*Phi';
    
    % k=k(x,x) of the corresponding stationary kernel.
    k = v;
    
    % normalize
    evs = k/kMax*e;
    
    if nargout > 1
        d = -4*(pi*F*rho*l).^2.*e;
        gevs = v*P*(d/sum(e)-e*sum(d)/(sum(e)^2));
    end
    
end