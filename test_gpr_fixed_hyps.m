
fname = 'Data/data_24-Oct-2017_with_gp.mat';
load(fname);

xstar = linspace(min(x),max(x),101)';

clear param;

param.domain = [-10,20];
param.m = 101;
param.sigma_eps = .1;
param.ev.l = 1;
param.ev.sigma = 1;

param.ef.fh = @pbc_ef;
param.ev.fh = @pbc_ev_normal;
param.ev.fh_L = @pbc_ev_normal_L;
% param.k.fh_L = @pse_L;


[mu,Sigma,L] = gpr_fixed_hyps(x,y,xstar,param);

fprintf('L =%f\n',L);

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
