function K = approx_covfunc(x1,x2,q,fh_eigenval,fh_eigenfunc,domain,kparam)

    evs = fh_eigenval(q,domain,kparam);
    Phi1 = fh_eigenfunc(x1,q,domain);
    Phi2 = fh_eigenfunc(x2,q,domain);
    
    K = Phi1*diag(evs)*Phi2';
end