function K = approx_covfunc(x1,x2,q,fh_eigenval,fh_eigenfunc,param)

    evs = fh_eigenval(q,param);
    Phi1 = fh_eigenfunc(x1,q,param);
    Phi2 = fh_eigenfunc(x2,q,param);
    
    K = Phi1*diag(evs)*Phi2';
end