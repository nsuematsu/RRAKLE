function evs = pse_eigenvalues(q,param)
% 周期的二乗指数型共分散関数の固有値列 [κ_0, κ_1, ..., κ_(q-1)] を返す関数．
%
% param: 共分散関数のパラメータを収めた構造体で，ここでは，
%   param.domain = [a,b]
%   kparam.l: スケールパラメータ，
%   kparam.v: 係数
%
% evs = [κ_0, κ_1, ..., κ_(q-1)]

    rho = 1/(param.domain(2)-param.domain(1));
    
    nu = floor((1:q)/2);
    x = 1/(2*pi*rho*param.l)^2;
    evs = param.v/rho*besseli(nu,x)/exp(x);
end