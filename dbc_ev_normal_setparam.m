function param = dbc_ev_normal_setparam(t,param)
% 周期条件のもとでの正規分布に従う固有値のパラメータ t を param に
% 反映させる関数．
%
%   θ_1 = log(σ_ε)
%   θ_2 = log(σ)
%   θ_3 = log(l)
% としているので以下の通り：
    param.sigma_eps = exp(t(1));
    param.ev.sigma = exp(t(2));
    param.ev.l = exp(t(3));
end