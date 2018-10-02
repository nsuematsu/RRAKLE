a = -10; b = 10;

clear param;

param.domain = [a,b];
param.m = 101;

param.ev.fh = @pbc_ev_studentT;
param.ev.sigma = 1;
param.ev.nu = 3;
param.ev.l = 1.5;

param.ef.fh = @pbc_ef;

x1 = linspace(a,b,101);
x2 = a;

[C,evs,Phi] = approx_covfunc(x1,x2,param);

fprintf('sum(evs)=%f and (b-a)=%f\n',sum(evs),b-a)

figure;
subplot(2,2,1)
param.k.sigma = param.ev.sigma;
param.k.nu = param.ev.nu;
param.k.l = param.ev.l;
Cmat = matern_k(x1,a,param);
plot(x1,C,x1,Cmat,'--')
xlabel('$x$','Interpreter','latex')
ylabel('$k(x,a)$','Interpreter','latex')
legend('Student-t','Matern');

x1b = linspace(a,a+0.1*(b-a),101);
Cb = approx_covfunc(x1b,x2,param);
Cbmat = matern_k(x1b,a,param);
subplot(2,2,2)
plot(x1b,Cb,x1b,Cbmat,'--')
xlabel('$x$','Interpreter','latex')
ylabel('$k(x,a)$','Interpreter','latex')
legend('Student-t','Matern')

subplot(2,2,3)
K = diag(sqrt(evs));
A = randn(param.m,5);
F = Phi*K*A;
plot(x1,F);
title('Samples')