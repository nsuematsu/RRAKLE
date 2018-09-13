function evs = discreteT_eigenvalues_for_periodic(param,k,b)

    a = param.domain(1); b = param.domain(2);
    m = param.m;
    k = param.ev.k;
    b = param.ev.b;
    
    P = b-a;
    
    J = floor(((0:m-1)+1)/2);
    Z = 8*hypergeom([1,-1i*b-k/2,1i*b-k/2],[1-1i*b+k/2,1+1i*b+k/2],1)...
        /((4*b^2+k^2)*pochhammer(1-k/2-1i*b,k)*pochhammer(1-k/2+1i*b,k))...
        -bare_t(0,k,-k/2,b);
    evs = real(P*bare_t(J,k,-k/2,b)/Z);
end

function v = bare_t(j,k,a,b)
    v = 1./((b^2 + (a+j).^2).*pochhammer(1+a-1i*b+j,k)...
        .*pochhammer(1+a+1i*b+j,k));
end