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
    
    param = param.fh_setparam(t);
    % param.sigma_eps = se;
    % param.ev.sigma = exp(t(2));
    % param.ev.l = exp(t(3));
    
    Phi = param.ef.fh(x,param);
    evs = param.ev.fh(param);
    Lam = diag(evs);
    
    % 不完全な reduced rank 近似（チェック用）
    % Sig = Phi*diag(evs)*Phi' + se^2*eye(n);
    % L = .5*(n*log(2*pi) + log(det(Sig)) + y'*(Sig\y));
    
    % reduced rank 近似
    A = Phi'*Phi*Lam+se^2*eye(m);
    L = .5*n*log(2*pi)+(n-m)*log(se)+.5*log(det(A));
    v = Phi'*y;
    L = L + .5/se^2*(y'*y-v'*(Lam/A)*v);
%     B = Phi'*Phi+se^2*diag(1./evs);
%     R = chol(B);    
%     v = R'\(Phi'*y);
%     L = L + + .5/se^2*(y'*y-v'*v);

    if nargout > 1
        
    end

end