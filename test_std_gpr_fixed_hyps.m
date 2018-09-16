
fname = 'Data/data_24-Oct-2017_with_gp.mat';
load(fname);

xstar = linspace(min(x),max(x),101)';

clear param;
param.sigma_eps = .1;
param.k.fh = @se_k;
param.k.fh_L = @se_L;
param.k.l = 1;
param.k.sigma = 1;

[mu,Sigma,L] = std_gpr_fixed_hyps(x,y,xstar,param);

fprintf('L=%f\n',L);

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
