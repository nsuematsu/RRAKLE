function evs = geo_eigenvalues_for_periodic(param,omega)
% 周期境界条件下で，固有値を幾何数列（指数関数）に基づいて設定するときの
% 固有値列を返す関す．
%  κ_j = A exp(-ω floor((j+1)/2))
% A は，
%   κ_0 + κ_1 + .... = ρ^(-1) = b - a
% を満たすように決める．
% omega (ω) は，正のパラメータで，大きいほど固有値は急激に小さくなる．
% つまり，より低周波成分よりのスペクトラムとなる．

    A = param.v*(param.domain(2)-param.domain(1))*(exp(omega)-1)/(exp(omega)+1);
    evs = A*exp(-omega*floor(((0:param.m-1)+1)/2));
end