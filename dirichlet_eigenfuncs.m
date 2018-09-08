function Phi = dirichlet_eigenfuncs_eval(x,q,param)

    a = param.domain(1); b = param.domain(2);
    rho = 1/(b-a);
    x = rho*(x(:)-a);
    n = length(x);
    
    aa = 1:q;
    bb = kron(x,aa);
    Phi = sqrt(2*rho)*sin(bb*pi);
end