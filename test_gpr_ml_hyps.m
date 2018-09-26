% run max_marginal_likelihood

fname = 'Data/data_24-Oct-2017_with_gp.mat';
load(fname);

xdag = linspace(min(x),max(x),101)';

t0 = [log(2),log(2),log(2)];
% t0 = [log(2),log(2)];

clear param;

param.domain = [-10,20];
param.m = 101;

param.ef.fh = @pbc_ef;

param.ev.fh = @pbc_ev_normal;
param.ev.fh_L = @pbc_ev_normal_L;

param.fh_setparam = @pbc_ev_normal_setparam;
opts.t0 = t0;
opts.use_gradient = false;

[mu,Sigma,t,L] = gpr_ml_hyps(x,y,xdag,param,opts);

fprintf('L=%f\n',L);

figure
std = sqrt(diag(Sigma));
upper = mu+3*std;
lower = mu-3*std;
area_x = [xdag(:)',fliplr(xdag(:)')];
area_y = [upper(:)',fliplr(lower(:)')];
fill(area_x,area_y,.75*ones(1,3))
hold on
plot(xdag,mu,x,y,'ro')
hold off

fprintf('Estimated σ_ε = %f\n', exp(t(1)));
fprintf('Estimated σ_k = %f\n', exp(t(2)));
fprintf('Estimated l = %f\n', exp(t(3)));
