% run max_marginal_likelihood

fname = 'Data/data_24-Oct-2017_with_gp.mat';
load(fname);

xstar = linspace(min(x),max(x),101)';



% use se
opts.t0 = [log(2),log(2),log(2)];
opts.use_gradient = true;
param.k.fh = @se_k;
param.k.fh_L = @se_L;
param.k.fh_setparam = @se_setparam;

% % use pse
% opts.t0 = [log(2),log(2),log(2)];
% opts.use_gradient = true;
% param.domain = [-10,20];
% param.k.fh = @pse_k;
% param.k.fh_L = @pse_L;
% param.k.fh_setparam = @pse_setparam;

% % use matern (opts.use_gradient must be false)
% opts.t0 = [log(2),log(2),log(2),log(10)];
% opts.use_gradient = false;
% param.k.fh = @matern_k;
% param.k.fh_L = @matern_L;
% param.k.fh_setparam = @matern_setparam;

[mu,Sigma,t] = std_gpr_ml_hyps(x,y,xstar,param,opts);

L = param.k.fh_L(t,x,y,param);
fprintf('L=%f\n',L)

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
if length(t)>3
    fprintf('Estimated nu = %f\n', exp(t(4)));
end
