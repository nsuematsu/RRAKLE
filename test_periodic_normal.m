a = 0; b = 10;

param.domain = [a,b];
param.m = 100;

param.ev.fh = @normal_eigenvalues_for_periodic;
param.ev.sigma = 1;
param.ev.l = 1;

param.ef.fh = @periodic_eigenfuncs;

x1 = linspace(a,b,101);
x2 = a;

[C,evs,Phi] = approx_covfunc(x1,x2,param);

fprintf('sum(evs)=%f and (b-a)=%f\n',sum(evs),b-a)

figure;
subplot(2,2,1)
param.k.sigma = 1;
param.k.l = param.ev.l;
Cpse = pse_k(x1,a,param);
Cse = se_k(x1,a,param);
plot(x1,C,x1,Cpse,x1,Cse,'--')
xlabel('$x$','Interpreter','latex')
ylabel('$k(x,a)$','Interpreter','latex')
legend('Normal','PSE','SE');

x1b = linspace(a,a+0.1*(b-a),101);
Cb = approx_covfunc(x1b,x2,param);
Cbpse = pse_k(x1b,a,param);
Cbse = pse_k(x1b,a,param);
subplot(2,2,2)
plot(x1b,Cb,x1b,Cbpse,x1b,Cbse,'--')
xlabel('$x$','Interpreter','latex')
ylabel('$k(x,a)$','Interpreter','latex')
legend('Normal','PSE','SE')

subplot(2,2,3)
K = diag(sqrt(evs));
A = randn(param.m,5);
F = Phi*K*A;
plot(x1,F);
title('Samples')