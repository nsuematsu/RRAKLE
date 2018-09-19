function [L,K,S] = pse_L(t,x,y,param)
%
% Inputs:
%   t(1) = θ_1 = log(σ_ε)
%   t(2) = θ_2 = log(σ_k)
%   t(3) = θ_3 = log(l)

    se = exp(t(1));
    sk = exp(t(2));
    l = exp(t(3));
    
    x = x(:); y = y(:);
    n = length(x);
    param.k.sigma = sk;
    param.k.l = l;
    K = pse_k(x,x,param);
    S = K + se^2*eye(n);
    L = 0.5*( n*log(2*pi)+log(det(S) )+y'*(S\y));
end