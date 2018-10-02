function [L,gradL] = pbc_ev_pse_L(t,x,y,param)
% 周期境界条件の下で，与えられたθでの
% 負の対数周辺尤度の reduced rank 近似を求める関数.
%
% Inputs:
%   t(1) = θ_1 = log(σ_ε)
%   t(2) = θ_2 = log(σ)
%   t(3) = θ_3 = log(l)

    m = param.m;
    n = length(x);
    se = exp(t(1));
    
    param = pbc_ev_pse_setparam(t,param);
    % param.sigma_eps = se;
    % param.ev.sigma = exp(t(2));
    % param.ev.l = exp(t(3));
    
    Phi = param.ef.fh(x,param);
    [evs,gevs] = param.ev.fh(param);
    Lam = diag(evs);
    
    % 不完全な reduced rank 近似（チェック用）
    % Sig = Phi*diag(evs)*Phi' + se^2*eye(n);
    % L = .5*(n*log(2*pi) + log(det(Sig)) + y'*(Sig\y));
    
    % reduced rank 近似
    A = Phi'*Phi*Lam+se^2*eye(m);
    Binv = Lam/A;
    C = A\(Phi'*Phi);
    v = Phi'*y;
    u = Binv*v;
    w = A\v;
    L = .5*n*log(2*pi)+(n-m)*log(se)+.5*log(det(A)) ...
        + .5/se^2*(y'*y-v'*u);

    if nargout > 1
        % dL/dθ1
        dL1dt1 = n-m;
    
        dL2dt1 = se^2*trace(inv(A));        
        
        dL3dt1 = -se^(-2)*(y'*y-v'*u) + u'*w;
        
        % dL/dθ2
        dL2dt2 = trace(C*Lam);
        
        dL3dt2 = -(v'*Binv)*w;        
        
        % dL/dθ3
        dLamdt3 = diag(gevs);
        
        dL2dt3 = 0.5*trace(C*dLamdt3);
        
        dL3dt3 = -0.5*w'*dLamdt3*w;
        
        % 
        gradL = [...
            dL1dt1 + dL2dt1 + dL3dt1; ...
            dL2dt2 + dL3dt2; ...
            dL2dt3 + dL3dt3 ...
            ];
        
    end

end