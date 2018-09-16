function evs = pbc_ev_studentT(param)
%
% Required param members:
%   param.m
%   param.ev.tau: Student-t's degree of freedom.
%   param.ev.scale: Student-t's scale parameter.

    a = param.domain(1); b = param.domain(2);
    m = param.m;
    
    P = b-a;
        
    v = param.ev.sigma^2;
    tau = param.ev.tau;
    sigma = param.ev.scale;
    
    S = floor(((0:m-1)+1)/2);
    
    % Make a tClocationScaleDistribution object.
    td = makedist('tLocationScale','mu',0,'sigma',sigma,'nu',tau);
    
    evs = pdf(td,S);

    % normalize
    evs = v*P*evs/sum(evs);
end