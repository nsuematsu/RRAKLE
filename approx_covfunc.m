function [K,evs,Phi1,Phi2] = approx_covfunc(x1,x2,q,fh_ev,fh_ef,param)

    evs = fh_ev(q,param);
    Phi1 = fh_ef(x1,q,param);
    Phi2 = fh_ef(x2,q,param);
    
    K = Phi1*diag(evs)*Phi2';
end