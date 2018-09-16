function evs = pbc_ev_pse(param)
% 周期的二乗指数型共分散関数の固有値列 [κ_0, κ_1, ..., κ_(m-1)] を返す関数．
%
% param: 共分散関数のパラメータを収めた構造体で，ここでは，
%   param.domain = [a,b]
%   param.l: スケールパラメータ，
%   param.v: 係数
%
% evs = [κ_0, κ_1, ..., κ_(q-1)]

    a = param.domain(1); b = param.domain(2);
    m = param.m;

    rho = 1/(b-a);
    
    nu = floor((1:m)/2);
    x = 1/(2*pi*rho*param.l)^2;
    evs = param.v/rho*besseli(nu,x)/exp(x);
end