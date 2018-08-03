function Phi = periodic_eigenfuncs_eval(x,q,param)
%
% x: n�~1�x�N�g��
% q: ��_0(x) ���� ��_{q-1}(x) �̒l�����߂�D
% param: �p�����[�^
%    param.domain: [a,b]
%
% Phi = [��_0(x(1)) ��_1(x(1)) ... ��_{q-1}(x(1));
%        ��_0(x(2)) ��_1(x(2)) ... ��_{q-1}(x(2));
%         .
%         .
%         .
%        ��_0(x(n)) ��_1(x(n)) ... ��_{q-1}(x(n))]

    a = param.domain(1); b = param.domain(2);
    rho = 1/(b-a);
    x = rho*(x(:)-a);
    n = length(x);
    Phi = ones(n,q)*sqrt(rho);
    
    aa = (1:q-1) + mod(1:q-1,2);
    bb = kron(x,aa);
    bb = bb + .5*repmat(mod(2:q,2),n,1);
    Phi(:,2:end) = sqrt(2*rho)*cos(bb*pi);
end