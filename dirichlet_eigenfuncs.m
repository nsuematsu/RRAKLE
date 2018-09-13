function Phi = dirichlet_eigenfuncs(x,param)

    a = param.domain(1); b = param.domain(2);
    m = param.m;
    
    rho = 1/(b-a);
    x = rho*(x(:)-a);
    
    aa = 1:m;
    bb = kron(x,aa);
    Phi = sqrt(2*rho)*sin(bb*pi);
end