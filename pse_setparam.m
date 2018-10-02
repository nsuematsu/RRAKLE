function param = pse_setparam(t,param)

    param.sigma_eps = exp(t(1));
    param.k.sigma = exp(t(2));
    param.k.l = exp(t(3));

end