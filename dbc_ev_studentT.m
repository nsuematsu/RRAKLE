function [evs,gevs] = dbc_ev_studentT(param)
%
% Required param members:
%   param.m
%   param.ev.tau: Student-t's degree of freedom.
%   param.ev.scale: Student-t's scale parameter.

    m = param.m;
    a = param.domain(1);
    b = param.domain(2);
    P = b-a;
    rho = 1/P;
    
    % 周期は 1/2, 1, 3/2, 2,... であることに注意
    F = .5*(1:m);
    
    v = param.ev.sigma^2;
    nu = param.ev.nu;
    l = param.ev.l;
      
    % Make a tLocationScaleDistribution object.
    tau = 2*nu; % DOF of t-distribution
    sigma = 1/(2*pi*rho*l); % scale of t-distribution
    td = makedist('tLocationScale','mu',0,'sigma',sigma,'nu',tau);
    
    e = pdf(td,F);

    % obatain max{k(x,x)}=k(c,c) when sum(κ_j)=1.
    c = mean(param.domain); % The center of the domain.
    Phi = param.ef.fh(c,param);
    kMax = Phi*diag(e)*Phi';
    
    % k(x,x) of the corresponding stationary kernel.
    k = v;
    
    % normalize
    evs = k/kMax*e;
    
    if nargout > 1
        a = 2*(pi*rho*l*F).^2;
        b = -sqrt(2)*pi*rho*l*nu^(nu+1)/beta(nu,.5);
        d1 = b*(nu+a).^(-nu-3/2).*(2*a-1);
        g1 = v*P*(d1/sum(e)-e*sum(d1)/(sum(e)^2));
        
        b = pi*rho*l*(nu^(nu+1))/(sqrt(2)*beta(nu,.5));
        d2 = b*(nu+a).^(-nu-3/2) ...
            .* (-1+2*a-2*(nu+a).*(log(nu+a)-log(nu)+psi(0,nu)-psi(0,nu+1/2)));
        g2 = v*P*(d2/sum(e)-e*sum(d2)/(sum(e)^2));
        
        gevs = [g1(:)'; g2(:)'];

    end
end