function [L,C,S] = se_L(t,x,y)
%
% Inputs:
%   t(1) = θ_1 = log(σ_ε)
%   t(2) = θ_2 = log(σ_k)
%   t(3) = θ_3 = log(l)

    se = exp(t(1));
    sk = exp(t(2));
    l = exp(t(3));
    
    x = x(:); y = y(:);
    param.k.sigma = sk;
    param.k.l = l;
    C = se_k(x,x,param);
    S = C + se^2*eye(size(x,1));
    L = 0.5*( size(x,1)*log(2*pi)+log(det(S) )+y'*(S\y));
end