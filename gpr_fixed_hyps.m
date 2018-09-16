function [mu,Sigma,L] = gpr_fixed_hyps(x,y,xstar,param)
% ガウス過程回帰を固定ハイパーパラメータで行う関数．
   
    m = param.m;
    se = param.sigma_eps;
    evs = param.ev.fh(param);
    Phi_x = param.ef.fh(x,param);
    Phi_xstar = param.ef.fh(xstar,param);
        
    K = diag(evs);
    Sigma_beta = se^2*K/((Phi_x'*Phi_x)*K+se^2*eye(m));
    mu_beta = Sigma_beta*(se^(-2)*Phi_x'*y);
    mu = Phi_xstar*mu_beta;
    Sigma = Phi_xstar*Sigma_beta*Phi_xstar';
    
    if nargout > 2
        t = [param.sigma_eps,param.ev.sigma,param.ev.l];
        L = param.k.fh_L(t,x,y);
    end
end