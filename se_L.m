function [L,gradL] = se_L(t,x,y,param)
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
    
    K = se_k(x,x,param);
    Sig = K + se^2*eye(n);
    L = 0.5*(n*log(2*pi)+log(det(Sig))+y'*(Sig\y));
    
    if nargout > 1
        gradL = zeros(3,1);
        Sig_rd_y = Sig\y;
    
        % dL/dt1
        dSig = 2*se^2*eye(length(x));
        gradL(1) = dLdt();
    
        % dL/dt2
        dSig = 2*Sig;
        gradL(2) = dLdt();
    
        % dL/dt3
        diff = repmat(x,1,length(x));
        S = (diff - diff').^2;    
        dSig = l^(-2)*S.*K;
        gradL(3) = dLdt();
    end

    function dL = dLdt()
        dL = 0.5*(trace(Sig\dSig)-Sig_rd_y'*dSig*Sig_rd_y);
    end
end