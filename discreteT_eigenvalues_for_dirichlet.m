function evs = discreteT_eigenvalues_for_dirichlet(param)
% 
    a = param.domain(1); b = param.domain(2);
    m = param.m;
  
    v = param.ev.sigma^2;
    k = param.ev.k;
    bb = param.ev.b;
        
    P = b-a;
    J = 1:m;

    % normalize so that its infinite sum to 1, first.
    Z = 4*hypergeom([1,-1i*bb-k/2+1,1i*bb-k/2+1],[2-1i*bb+k/2,2+1i*bb+k/2],1)...
        /((4+4*bb^2-4*k+k^2)*pochhammer(2-k/2-1i*bb,k)*pochhammer(2-k/2+1i*bb,k));
    evs = real(bare_t(J,k,-k/2,bb)/Z);
    
    % obatain max{k(x,x)}=k(c,c) when sum(κ_j)=1.
    c = .5*(a+b);
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