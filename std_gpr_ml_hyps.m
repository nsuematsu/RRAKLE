function [mu,Sigma,t] = std_gpr_ml_hyps(x,y,xstar,param,opts)
% GPR with maximum marginal likelihood hyperparameters.

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
    % fh = @(t)L_gradL(t,x,y,param);
    t = fminunc(@L_gradL,t0,optimOpts,x,y,param);
    
    % Calc mu and Sigma
    param = param.fh_setparam(t,param);
    se = param.sigma_eps; % σ_ε    
    Kxx = param.k.fh(x,x,param);
    S = Kxx+se^2*eye(size(x,1));
    Kxxdag = param.k.fh(x,xstar,param);
    Kxdagxdag = param.k.fh(xstar,xstar,param);
    tmp = Kxxdag'/S;
    
    mu = tmp*y;
    Sigma = Kxdagxdag - tmp*Kxxdag;    
end

function [L,gradL] = L_gradL(t,x,y,param)
    [L,C,S] = param.k.fh_L(t,x,y);
    if nargout > 1
        gradL = param.k.fh_gradL(t,x,y,C,S);
    end
end
