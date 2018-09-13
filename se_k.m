function K = se_k(x1,x2,param)
% 二乗指数型共分散関数

    l = param.k.l;
    v = param.k.sigma^2;

    n1 = numel(x1);
    n2 = numel(x2);
    
    X1 = repmat(x1(:),1,n2);
    X2 = repmat(x2(:),1,n1);
    d2 = (X1 - X2').^2;
    K = v*exp(-1/(2*l^2)*d2);
    
end