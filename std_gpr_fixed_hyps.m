function [mu,Sigma,nlml] = std_gpr_fixed_hyps(x,y,xstar,param)
% 標準的ガウス過程回帰を固定ハイパーパラメータで行う関数．

    x = x(:); y = y(:);
    n = length(x);

    K_x_x = param.k.fh(x,x,param);
    K_xstar_x = param.k.fh(xstar,x,param);
    K_xstar_xstar = param.k.fh(xstar,xstar,param);
    
    R = chol(K_x_x + param.sigma_eps.^2 * eye(n));
    A = K_xstar_x/R/R';
    mu = A*y;
    Sigma = K_xstar_xstar-A*K_xstar_x';
    
    if nargout > 2 
        % Given θ, calc negative log marginal likelihood
        t = log([param.sigma_eps,param.k.sigma,param.k.l]); % θ
        nlml = param.k.fh_L(t,x,y);
    end

end