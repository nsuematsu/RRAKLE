function param = pbc_ev_pse_setparam(t,param)

    param.sigma_eps = exp(t(1));
    param.ev.sigma = exp(t(2));
    param.ev.l = exp(t(3));
end