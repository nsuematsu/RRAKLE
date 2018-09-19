function [K,evs,Phi1,Phi2] = approx_covfunc(x1,x2,param)

    evs = param.ev.fh(param);
    Phi1 = param.ef.fh(x1,param);
    Phi2 = param.ef.fh(x2,param);
    
    K = Phi1*diag(evs)*Phi2';
end