function evs = dbc_ev_studentT(param)
%
% Required param members:
%   param.m
%   param.ev.tau: Student-t's degree of freedom.
%   param.ev.scale: Student-t's scale parameter.

    m = param.m;
  
    v = param.ev.sigma^2;
    tau = param.ev.tau;
    sigma = param.ev.scale;
    J = 1:m;
    
    % Make a tClocationScaleDistribution object.
    td = makedist('tLocationScale','mu',0,'sigma',sigma,'nu',tau);
    
    evs = pdf(td,J);

    % obatain max{k(x,x)}=k(c,c) when sum(κ_j)=1.
    c = mean(param.domain); % The center of the domain.
    Phi = param.ef.fh(c,param);
    kMax = Phi*diag(evs)*Phi';
    
    % k(x,x) of the corresponding stationary kernel.
    k = v;
    
    % normalize
    evs = k/kMax*evs;
end

function v = bare_t(j,k,a,b)
    v = 1./((b^2 + (a+j).^2).*pochhammer(1+a-1i*b+j,k)...
        .*pochhammer(1+a+1i*b+j,k));
end