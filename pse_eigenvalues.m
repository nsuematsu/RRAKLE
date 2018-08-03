function evs = pse_eigenvalues(q,param)
% �����I���w���^�����U�֐��̌ŗL�l�� [��_0, ��_1, ..., ��_(q-1)] ��Ԃ��֐��D
%
% param: �����U�֐��̃p�����[�^�����߂��\���̂ŁC�����ł́C
%   param.domain = [a,b]
%   kparam.l: �X�P�[���p�����[�^�C
%   kparam.v: �W��
%
% evs = [��_0, ��_1, ..., ��_(q-1)]

    rho = 1/(param.domain(2)-param.domain(1));
    
    nu = floor((1:q)/2);
    x = 1/(2*pi*rho*param.l)^2;
    evs = param.v/rho*besseli(nu,x)/exp(x);
end