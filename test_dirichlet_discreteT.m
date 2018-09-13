a = 0; b = 10;
c = 0.5*(a+b);
m = 200;

clear param;

param.domain = [a,b];
param.m = m;

param.ev.fh = @discreteT_eigenvalues_for_dirichlet;
param.ev.sigma = 1;
param.ev.k = 3/2;
param.ev.b = 10;

param.ef.fh = @dirichlet_eigenfuncs;

x1 = linspace(a,b,101);

[C,evs,Phi] = approx_covfunc(x1,x1,param);

fprintf('sum(evs)=%f and (b-a)=%f\n',sum(evs),b-a)

figure;
subplot(2,2,1)
param.k.sigma = param.ev.sigma;
param.k.nu = 2*param.ev.k+1;
param.k.l =1/sqrt(param.ev.k^4/4+param.ev.b^2);
Cmat = matern_k(x1,c,param);
index_c = find(x1==c);
plot(x1,C(index_c,:),x1,Cmat,'--');
% plot(x1,C(index_c,:));
xlabel('$x$','Interpreter','latex')
ylabel('$k(x,c)$ or $k_\mathrm{Matern}(x,c)$','Interpreter','latex')
% ylabel('$k(x,c)$','Interpreter','latex')
subplot(2,2,2)
[X,Y]=meshgrid(x1);
mesh(X,Y,C);

x1b = linspace(c-0.1*(b-a),c+0.1*(b-a),101);
Cb = approx_covfunc(x1b,c,param);
subplot(2,2,3)
plot(x1b,Cb)
xlabel('$x$','Interpreter','latex')
ylabel('$k(x,c)$','Interpreter','latex')

subplot(2,2,4)
K = diag(sqrt(evs));
A = randn(param.m,5);
F = Phi*K*A;
plot(x1,F);
title('Samples')