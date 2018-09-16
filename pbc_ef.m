function Phi = pbc_ef(x,param)
% pbc_ef returns the values of the m eigenfunctions for the Periodic BC.
%
% x: n×1ベクトル
% q: φ_0(x) から φ_{q-1}(x) の値を求める．
% param: パラメータ
%    param.domain: [a,b]
%
% Phi = [φ_0(x(1)) φ_1(x(1)) ... φ_{q-1}(x(1));
%        φ_0(x(2)) φ_1(x(2)) ... φ_{q-1}(x(2));
%         .
%         .
%         .
%        φ_0(x(n)) φ_1(x(n)) ... φ_{q-1}(x(n))]

    a = param.domain(1); b = param.domain(2);
    m = param.m;
    
    rho = 1/(b-a);
    x = rho*(x(:)-a);
    n = length(x);
    Phi = ones(n,m)*sqrt(rho);
    
    aa = (1:m-1) + mod(1:m-1,2);
    bb = kron(x,aa);
    bb = bb + .5*repmat(mod(2:m,2),n,1);
    Phi(:,2:end) = sqrt(2*rho)*cos(bb*pi);
end