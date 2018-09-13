function [mu,Sigma] = std_gpr_fixed_hyps(x,y,xstar,param)
% 標準的ガウス過程回帰を固定ハイパーパラメータで行う関数．

    n = length(x);

    k_x_x = param.fh_k(x,x,param);
    k_xstar_x = param.fh_k(xstar,x,param);
    
    tmp = k_xstar_x / (k_x_x + param.s_eps.^2 * eye(n));
    mu = tmp * y(:);
    Sigma = tmp * k_xstar_x';

end