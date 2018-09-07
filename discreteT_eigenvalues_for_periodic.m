function evs = discreteT_eigenvalues_for_periodic(q,param,k,b)
    L = param.domain(2) - param.domain(1);
    J = floor(((0:q-1)+1)/2);
    Z = 8*hypergeom([1,-1i*b-k/2,1i*b-k/2],[1-1i*b+k/2,1+1i*b+k/2],1)...
        /((4*b^2+k^2)*pochhammer(1-k/2-1i*b,k)*pochhammer(1-k/2+1i*b,k))...
        -bare_t(0,k,-k/2,b);
    evs = real(L*bare_t(J,k,-k/2,b)/Z);
end

function v = bare_t(j,k,a,b)
    v = 1./((b^2 + (a+j).^2).*pochhammer(1+a-1i*b+j,k)...
        .*pochhammer(1+a+1i*b+j,k));
end