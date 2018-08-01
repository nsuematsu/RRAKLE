function evs = pse_eigenvalues(q,domain,kparam)
% 周期的二乗指数型共分散関数の固有値列 [κ_0, κ_1, ..., κ_(q-1)] を返す関数．
%
% domain: [a,b] もしくは ρ=1/(b-a)
% kparam: 共分散関数のパラメータを収めた構造体で，ここでは，
%   kparam.l: スケールパラメータ，
%   kparam.v: 係数
%
% evs = [κ_0, κ_1, ..., κ_(q-1)]

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