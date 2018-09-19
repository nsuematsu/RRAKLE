param.m = 101;
param.domain = [-10,10];

param.k.fh = @pse_k;
param.k.sigma = 1;
param.k.l = 1;

param.ef.fh = @pbc_ef;
param.ev.fh = @pbc_ev_pse;
param.ev.sigma = 1;
param.ev.l = 1;

x = linspace(-5,5,21);
Kxx = param.k.fh(x,x,param);

Kxx2 = approx_covfunc(x,x,param);

figure;
subplot(1,2,1)
mesh(Kxx)
subplot(1,2,2)
mesh(Kxx2)