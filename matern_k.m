function C = matern_k(x1,x2,param)
% Matern convlution function

    v = param.k.sigma^2;
    nu = param.k.nu;
    l = param.k.l;

    n1 = numel(x1);
    n2 = numel(x2);
    
    X1 = repmat(x1(:),1,n2);
    X2 = repmat(x2(:),1,n1);
    A = sqrt(2*nu)/l*abs(X1 - X2');
    C = 2^(1-nu)/gamma(nu)*(A).^nu;
    C = v*C.*besselk(nu,A);    
end