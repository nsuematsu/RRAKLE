function [evs,gevs] = pbc_ev_studentT(param)
%
% Required param members:
%   param.m
%   param.ev.nu: DOF of the corresponding Matern kernel.
%   param.ev.l: scale parameter of the corresponding Matern kernel.

    a = param.domain(1); b = param.domain(2);
    m = param.m;
    
    P = b-a;
    rho = 1/P;
        
    v = param.ev.sigma^2;
    nu = param.ev.nu;
    l = param.ev.l;
    
    S = floor(((0:m-1)+1)/2);
    
    % Make a tClocationScaleDistribution object.
    tau = 2*nu; % DOF of t-distribution
    sigma = 1/(2*pi*rho*l); % scale of t-distribution
    td = makedist('tLocationScale','mu',0,'sigma',sigma,'nu',tau);
    
    e = pdf(td,S);

    % normalize
    evs = v*P*e/sum(e);
    
    if nargout > 1        
        a = 2*(pi*rho*l*S).^2;
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