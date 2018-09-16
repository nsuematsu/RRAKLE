function [mu,Sigma,t] = std_gpr_ml_hyps(x,y,xstar,param,opts)
% GPR with maximum marginal likelihood hyperparameters.

    t0 = opts.t0; % Initial values
    
    % Options for the optimizer.
    optimOpts = optimoptions('fminunc',...
        'Algorithm','trust-region',...
        'SpecifyObjectiveGradient',true);
    
    % Perform the optimization.
    % fh = @(t)L_gradL(t,x,y,param);
    t = fminunc(@L_gradL,t0,optimOpts,x,y,param);
    
    % Calc mu and Sigma
    se = exp(t(1)); % σ_ε    
    param.k.sigma = exp(t(2));
    param.k.l = exp(t(3));
    C_x_x = param.k.fh(x,x,param);
    S = C_x_x+se^2*eye(size(x,1));
    C_x_xs = param.k.fh(x,xstar,param);
    C_xs_xs = param.k.fh(xstar,xstar,param);
    tmp = C_x_xs'/S;
    
    mu = tmp*y;
    Sigma = C_xs_xs - tmp*C_x_xs;    
end

function [L,gradL] = L_gradL(t,x,y,param)
    [L,C,S] = param.k.fh_L(t,x,y);
    if nargout > 1
        gradL = param.k.fh_gradL(t,x,y,C,S);
    end
end
