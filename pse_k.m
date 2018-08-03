function K = pse_k(x1,x2,param)

    rho = 1/(param.domain(2) - param.domain(1));
    n1 = numel(x1);
    n2 = numel(x2);
    
    X1 = repmat(x1(:),1,n2);
    X2 = repmat(x2(:),1,n1);
    S2 = sin(pi*rho*(X1 - X2')).^2;
    K = param.v*exp(-2/(2*pi*rho*param.l)^2*S2);
    
end