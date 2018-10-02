function [mu_dag,Sig_dag,t,L,output] = gpr_ml_hyps(x,y,xdag,param,opts)
% 提案手法によるガウス過程回帰を周辺尤度最大化で求めたハイパーパラメータで行う関数．

    t0 = opts.t0; % Initial values
        
    % Options for the optimizer.
    if opts.use_gradient
        optimOpts = optimoptions('fminunc',...
            'Algorithm','trust-region',...
            'SpecifyObjectiveGradient',true);
    else
        optimOpts = optimoptions('fminunc',...
            'Algorithm','quasi-newton',...
            'SpecifyObjectiveGradient',false);
    end
    
    % Perform the optimization.
    [t,L,~,output] = fminunc(param.ev.fh_L,t0,optimOpts,x,y,param);
    
    % Calc mu and Sigma
    param = param.ev.fh_setparam(t,param);
    
    m = param.m;
    se = param.sigma_eps;
    evs = param.ev.fh(param);
    Phi_x = param.ef.fh(x,param);
    Phi_xdag = param.ef.fh(xdag,param);

    Lam = diag(evs); % Λ
    Sig = se^2*Lam/((Phi_x'*Phi_x)*Lam+se^2*eye(m)); % Σ
    mu = Sig*(se^(-2)*Phi_x'*y); % μ
    mu_dag = Phi_xdag*mu;
    Sig_dag = Phi_xdag*Sig*Phi_xdag';    
end