% run max_marginal_likelihood

fname = 'Data/data_24-Oct-2017_with_gp.mat';
load(fname);

xstar = linspace(min(x),max(x),101)';

t0 = [log(2),log(2),log(2)];

param.k.fh = @se_k;
param.k.fh_L = @se_L;
param.k.fh_gradL = @se_gradL;
opts.t0 = t0;

[mu,Sigma,t] = std_gpr_ml_hyps(x,y,xstar,param,opts);

figure
std = sqrt(diag(Sigma));
upper = mu+3*std;
lower = mu-3*std;
area_x = [xstar(:)',fliplr(xstar(:)')];
area_y = [upper(:)',fliplr(lower(:)')];
fill(area_x,area_y,.75*ones(1,3))
hold on
plot(xstar,mu,x,y,'ro')
hold off

fprintf('Estimated σ_ε = %f\n', exp(t(1)));
fprintf('Estimated σ_k = %f\n', exp(t(2)));
fprintf('Estimated l = %f\n', exp(t(3)));
