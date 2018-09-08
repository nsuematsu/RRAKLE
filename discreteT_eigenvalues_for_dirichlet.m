function evs = discreteT_eigenvalues_for_dirichlet(q,param,k,b)
    L = param.domain(2) - param.domain(1);
    J = 1:q;
    Z = 4*hypergeom([1,-1i*b-k/2+1,1i*b-k/2+1],[2-1i*b+k/2,2+1i*b+k/2],1)...
        /((4+4*b^2-4*k+k^2)*pochhammer(2-k/2-1i*b,k)*pochhammer(2-k/2+1i*b,k));
    evs = real(L*bare_t(J,k,-k/2,b)/Z);
end

function v = bare_t(j,k,a,b)
    v = 1./((b^2 + (a+j).^2).*pochhammer(1+a-1i*b+j,k)...
        .*pochhammer(1+a+1i*b+j,k));
end