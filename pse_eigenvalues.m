function evs = pse_eigenvalues(q,domain,kparam)
% �����I���w���^�����U�֐��̌ŗL�l�� [��_0, ��_1, ..., ��_(q-1)] ��Ԃ��֐��D
%
% domain: [a,b] �������� ��=1/(b-a)
% kparam: �����U�֐��̃p�����[�^�����߂��\���̂ŁC�����ł́C
%   kparam.l: �X�P�[���p�����[�^�C
%   kparam.v: �W��
%
% evs = [��_0, ��_1, ..., ��_(q-1)]

    if numel(domain) == 2
        rho = 1/(domain(2)-domain(1));
    elseif numel(domain) == 1
        rho = domain;
    else
        error('domain must be [a,b] or 1/(b-a).')
    end
    
    nu = floor((1:q)/2);
    x = 1/(2*pi*rho*kparam.l)^2;
    evs = kparam.v/rho*besseli(nu,x)/exp(x);
end