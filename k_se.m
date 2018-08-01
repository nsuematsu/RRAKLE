function K = k_se(x1,x2,kparam)

    n1 = numel(x1);
    n2 = numel(x2);
    
    X1 = repmat(x1(:),1,n2);
    X2 = repmat(x2(:),1,n1);
    d2 = (X1 - X2').^2;
    K = exp(-1/(2*kparam.l^2)*d2);
    
end