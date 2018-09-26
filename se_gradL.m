function gradL = se_gradL(t,x,y,C,Sig)
% gradL_se returns a gradient vector of the negative log likelihood
% for SE kernel.
% Note that m(x)=0 is assumed in this implementation.
%
%
% Inputs:
%   t(1) = θ_1 = log(σ_ε)
%   t(2) = θ_2 = log(σ_k)
%   t(3) = θ_3 = log(l)

    se = exp(t(1));
    sk = exp(t(2));
    l = exp(t(3));
    
    x = x(:); y = y(:);

    if nargin < 5
        param.k.sigma = sk;
        param.k.l = l;
        C = se_k(x,x,param); % C_{x,x}
        Sig = C + exp(2*se)*eye(length(x)); % Σ
    end
    
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
    dSig = l^(-2)*S.*C;
    gradL(3) = dLdt();

    function dL = dLdt()
        dL = 0.5*(trace(Sig\dSig)-Sig_rd_y'*dSig*Sig_rd_y);
    end
end