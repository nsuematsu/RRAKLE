function [evs,gevs] = pbc_ev_pse(param)
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
    
    l = param.ev.l;
    v = param.ev.sigma^2;
    P = b-a;
    rho = 1/P;
    
    S = floor((1:m)/2);
    x = (2*pi*rho*l)^(-2);
    evs = v*P*besseli(S,x)/exp(x);
    
    if nargout > 1
        gevs = .5*v*P*x*exp(-x)* ...
               (besseli(S,x).*(1+S/x)-besseli(S-1,x));
        
    end
    
end