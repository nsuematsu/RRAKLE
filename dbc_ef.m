function Phi = dbc_ef(x,param)
% dbc_ef returns the values of the m eigenfunctions for the Dirichlet BC.
% x: n×1ベクトル
% Required param's members:
%    prram.m 
%    param.domain: [a,b]
%
% Phi = [φ_0(x(1)) φ_1(x(1)) ... φ_{m-1}(x(1));
%        φ_0(x(2)) φ_1(x(2)) ... φ_{m-1}(x(2));
%         .
%         .
%         .
%        φ_0(x(n)) φ_1(x(n)) ... φ_{m-1}(x(n))]

    a = param.domain(1); b = param.domain(2);
    m = param.m;
    
    rho = 1/(b-a);
    x = rho*(x(:)-a);
    
    aa = 1:m;
    bb = kron(x,aa);
    Phi = sqrt(2*rho)*sin(bb*pi);
end