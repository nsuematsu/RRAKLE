function evs = geo_eigenvalues_for_dirichlet(param,omega)
% ディリクレ境界条件下で，固有値を幾何数列（指数関数）に基づいて設定するときの
% 固有値列を返す関す．
%  κ_j = A exp(-ω (j-1))
% A は，
%   κ_1 + κ_2 + .... = ρ^(-1) = b - a
% を満たすように決める．
% omega (ω) は，正のパラメータで，大きいほど固有値は急激に小さくなる．
% つまり，より低周波成分よりのスペクトラムとなる．

    A = param.v*(param.domain(2)-param.domain(1))*(1-exp(-omega));
    evs = A*exp(-omega*(1:param.m));
end