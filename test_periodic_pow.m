a = -10; b = 10;

param.domain = [a,b];
param.v = 1;

x1 = linspace(a,b,101);
x2 = a;

fh_eigv = @(q,param) pow_eigenvalues_for_periodic(q,param,3);
fh_eigf = @periodic_eigenfuncs;

q = 1000;
K = approx_covfunc(x1,x2,q,fh_eigv,fh_eigf,param);

plot(x1,K)