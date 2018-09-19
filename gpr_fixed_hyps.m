function [mu_dag,Sigma_dag,L] = gpr_fixed_hyps(x,y,xdag,param)
% ガウス過程回帰を固定ハイパーパラメータで行う関数．
   
    m = param.m;
    se = param.sigma_eps;
    evs = param.ev.fh(param);
    Phi_x = param.ef.fh(x,param);
    Phi_xdag = param.ef.fh(xdag,param);

    Lam = diag(evs); % Λ
    Sig = se^2*Lam/((Phi_x'*Phi_x)*Lam+se^2*eye(m)); % Σ
    mu = Sig*(se^(-2)*Phi_x'*y); % μ
    mu_dag = Phi_xdag*mu;
    Sigma_dag = Phi_xdag*Sig*Phi_xdag';
    
    if nargout > 2
        t = log([param.sigma_eps,param.ev.sigma,param.ev.l]);
        L = param.ev.fh_L(t,x,y,param);
    end
end