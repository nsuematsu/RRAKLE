function evs = pow_eigenvalues_for_periodic(q,param,p)
% 周期境界条件下で，固有値をべき関数に基づいて設定するときの
% 固有値列を返す関す．
%  κ_j = A floor((j+3)/2)^{-p}
% A は，
%   κ_0 + κ_1 + .... = ρ^(-1) = b - a
% を満たすように決める．
% p は，1より大きなパラメータで，大きいほど固有値は急激に小さくなる．
% つまり，より低周波成分よりのスペクトラムとなる．

    A = param.v*(param.domain(2)-param.domain(1))/(2*zeta(p)-1);
    evs = A*(floor(((0:q-1)+3)/2)).^(-p);
end