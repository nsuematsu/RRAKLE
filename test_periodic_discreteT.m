a = 0; b = 10;

param.domain = [a,b];
param.m = 100;

param.ev.fh = @discreteT_eigenvalues_for_periodic;
param.ev.sigma = 1;
param.ev.k = 3/2; % ポッホハマー記号で計算しているので整数でなくてよい
param.ev.b = 4;

param.ef.fh = @periodic_eigenfuncs;


x1 = linspace(a,b,101);
x2 = a;

[K,evs,Phi] = approx_covfunc(x1,x2,param);

fprintf('sum(evs)=%f and (b-a)=%f\n',sum(evs),b-a)

subplot(2,2,1)
plot(x1,K)
xlabel('$x$','Interpreter','latex')
ylabel('$k(x,a)$','Interpreter','latex')

x1b = linspace(a,a+0.1*(b-a),101);
Kb = approx_covfunc(x1b,x2,param);
subplot(2,2,2)
plot(x1b,Kb)
xlabel('$x$','Interpreter','latex')
ylabel('$k(x,a)$','Interpreter','latex')

subplot(2,2,3)
K = diag(sqrt(evs));
A = randn(param.m,5);
F = Phi*K*A;
plot(x1,F);
title('Samples')