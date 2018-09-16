function [mu,Sigma,L] = std_gpr_fixed_hyps(x,y,xstar,param)
% 標準的ガウス過程回帰を固定ハイパーパラメータで行う関数．

    n = size(x,1);

    C_x_x = param.k.fh(x,x,param);
    C_xstar_x = param.k.fh(xstar,x,param);
    C_xstar_xstar = param.k.fh(xstar,xstar,param);
    
    tmp = C_xstar_x / (C_x_x + param.sigma_eps.^2 * eye(n));
    mu = tmp * y(:);
    Sigma = C_xstar_xstar - tmp*C_xstar_x';
    
    if nargout > 2
        t = [param.sigma_eps,param.k.sigma,param.k.l];
        L = param.k.fh_L(t,x,y);
    end

end