% run gpr_ml_hyps

% fname = 'Data/data_24-Oct-2017_with_gp.mat';
fname = 'data.mat';
load(fname);

xdag = linspace(min(x),max(x),101)';

clear param;

param.domain = [0,10];
param.m = 25;

param.ef.fh = @pbc_ef;

% type = 'pse';
% type = 'pbc_normal';
% type = 'pbc_studentT';
type = 'dbc_normal';
% type = 'dbc_studentT';
switch type
    case 'pse' % use pse
        opts.t0 = [log(2),log(2),log(2)];
        opts.use_gradient = true;
        param.ev.fh = @pbc_ev_pse;
        param.ev.fh_L = @pbc_ev_pse_L;
        param.ev.fh_setparam = @pbc_ev_pse_setparam;
    case 'pbc_normal' % use pbc_normal
        opts.t0 = [log(2),log(2),log(2)];
        opts.use_gradient = true;
        param.ev.fh = @pbc_ev_normal;
        param.ev.fh_L = @pbc_ev_normal_L;
        param.ev.fh_setparam = @pbc_ev_normal_setparam;
    case 'pbc_studentT' % use pbc_studentT
        opts.t0 = [log(2),log(2),log(2),log(10)];
        opts.use_gradient = true;
        param.ev.fh = @pbc_ev_studentT;
        param.ev.fh_L = @pbc_ev_studentT_L;
        param.ev.fh_setparam = @pbc_ev_studentT_setparam;
    case 'dbc_normal' % use dbc_normal
        opts.t0 = [log(2),log(2),log(2)];
        opts.use_gradient = true;
        param.ev.fh = @dbc_ev_normal;
        param.ev.fh_L = @dbc_ev_normal_L;
        param.ev.fh_setparam = @dbc_ev_normal_setparam;
    case 'dbc_studentT' % use dbc_StudentT
        opts.t0 = [log(2),log(2),log(2),log(10)];
        opts.use_gradient = true;
        param.ev.fh = @dbc_ev_studentT;
        param.ev.fh_L = @dbc_ev_studentT_L;
        param.ev.fh_setparam = @dbc_ev_studentT_setparam;
    otherwise
        error('Unknown type:%s',type)
end

[mu,Sigma,t,L,output] = gpr_ml_hyps(x,y,xdag,param,opts);

fprintf('L=%f\n',L);
fprintf('iterations=%d\n',output.iterations);

figure
std = sqrt(diag(Sigma));
upper = mu+3*std;
lower = mu-3*std;
area_x = [xdag(:)',fliplr(xdag(:)')];
area_y = [upper(:)',fliplr(lower(:)')];
fill(area_x,area_y,.75*ones(1,3))
hold on
plot(xdag,mu,x,y,'ro')
title(type,'Interpreter','none')
hold off

fprintf('Estimated σ_ε = %f\n', exp(t(1)));
fprintf('Estimated σ_k = %f\n', exp(t(2)));
fprintf('Estimated l = %f\n', exp(t(3)));
if length(t)>3
    fprintf('Estimated nu = %f\n', exp(t(4)));
end