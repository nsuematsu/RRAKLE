a = -10; b = 9;
c = 0.5*(a+b);
rho = 1/(b-a);
m = 200;

clear param;

param.domain = [a,b];
param.m = m;

param.ef.fh = @dbc_ef;

param.ev.fh = @dbc_ev_studentT;
param.ev.sigma = 1;
param.ev.nu = 3;
param.ev.l = 2;

x1 = linspace(a,b,101);

[C,evs,Phi] = approx_covfunc(x1,x1,param);

fprintf('sum(evs)=%f and (b-a)=%f\n',sum(evs),b-a)

figure;
subplot(2,2,1)
param.k.sigma = param.ev.sigma;
param.k.nu = param.ev.nu;
param.k.l = param.ev.l;
Cmat = matern_k(x1,c,param);
index_c = find(x1==c);
plot(x1,C(index_c,:),x1,Cmat,'--');
% plot(x1,C(index_c,:));
xlabel('$x$','Interpreter','latex')
ylabel('$k(x,c)$ and $k_\mathrm{Matern}(x,c)$','Interpreter','latex')
% ylabel('$k(x,c)$','Interpreter','latex')
xlim([x1(1),x1(end)])
subplot(2,2,2)
[X,Y]=meshgrid(x1);
mesh(X,Y,C);

x1b = linspace(c-0.1*(b-a),c+0.1*(b-a),101);
Cb = approx_covfunc(x1b,c,param);
Cbmat = matern_k(x1b,c,param);
subplot(2,2,3)
plot(x1b,Cb,x1b,Cbmat,'--')
xlabel('$x$','Interpreter','latex')
ylabel('$k(x,c)$ and $k_\mathrm{Matern}(x,c)$','Interpreter','latex')
xlim(x1b([1,end]))

subplot(2,2,4)
K = diag(sqrt(evs));
A = randn(param.m,5);
F = Phi*K*A;
plot(x1,F);
xlim([x1(1),x1(end)])

title('Samples')