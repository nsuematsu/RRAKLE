function [L,gradL] = pbc_ev_normal_L(t,x,y,param)
% 周期境界条件の下で，与えられたθでの
% 負の対数周辺尤度の reduced rank 近似を求める関数.
%
% 出力変数が2の時，Lのθに関する勾配も返す．この勾配は，
% ハイパーパラメータの学習に利用．
%
% Inputs:
%   t(1) = θ_1 = log(σ_ε)
%   t(2) = θ_2 = log(σ)
%   t(3) = θ_3 = log(l)

    m = param.m;
    n = length(x);
    se = exp(t(1));
    
    param = param.fh_setparam(t,param);
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
    v = Phi'*y;
    u = Binv*v;
    L = .5*n*log(2*pi)+(n-m)*log(se)+.5*log(det(A))+...
            .5/se^2*(y'*y-v'*u);
    
    % dL/dθ is required.
    if nargout > 1
        
        % dL/dθ1
        dL1dt1 = n-m;
    
        dL2dt1 = se^2*trace(inv(A));        
        
        % dL3dt1 = -se^(-2)*(y'*y-v'*u) + u'*(Lam\u);
        dL3dt1 = -se^(-2)*(y'*y-v'*u) + u'*(A\v);
        
        % dL/dθ2
        dL2dt2 = trace(A\(Phi'*Phi)*Lam);
        
        dL3dt2 = -(v'*Binv)*(A\v);        
        
        % dL/dθ3
        dLamdt3 = diag(gevs);
        
        dL2dt3 = trace(A\(Phi'*Phi)*dLamdt3);
        
        w = A\v;
        dL3dt3 = -0.5*w'*dLamdt3*w;
        
        % 
        gradL = [...
            dL1dt1 + dL2dt1 + dL3dt1; ...
            dL2dt2 + dL3dt2; ...
            dL2dt3 + dL3dt3 ...
            ];
    
    end 
end