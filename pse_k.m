function K = pse_k(x1,x2,param)
% 周期的二乗指数型共分散関数

    a = param.domain(1); b = param.domain(2);
    l = param.k.l;
    v = param.k.sigma^2;

    rho = 1/(b - a);
    n1 = numel(x1);
    n2 = numel(x2);
    
    X1 = repmat(x1(:),1,n2);
    X2 = repmat(x2(:),1,n1);
    S2 = sin(pi*rho*(X1 - X2')).^2;
    K = v*exp(-2/(2*pi*rho*l)^2*S2);
    
end