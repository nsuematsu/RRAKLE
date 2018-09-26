function [evs,gevs] = pbc_ev_normal(param)
% pbc_ev_normal は周期条件のもとでの正規分布に基づく固有値を返す．
% 
% 出力変数が2の時，固有値のパラメータに関する勾配も返す．
% 勾配と言っても，今パラメータが一つしかないので微分を返すだけ．
% より詳しくは θ=log(l) として，θに関する微分を返す．
% 
    a = param.domain(1); b = param.domain(2);
    m = param.m;
    v = param.ev.sigma^2;
    l = param.ev.l;
    P = b-a;
    rho = 1/P;

    S = floor(((0:m-1)+1)/2);
    tmp = exp(-2*(pi*S*rho*l).^2);
    
    % 無限項まで考慮した正規化
    % evs = v*tmp;
    % evs = P/theta3(exp(-2*(pi*rho*l).^2))*evs;
    
    % m項までで正規化
    evs = v*P/sum(tmp)*tmp;
    
    if nargout > 1 % 勾配も必要（勾配と言ってもここでは単なる部分）
        % とりあえず，正規化を無視した微分
        gevs = -4*(pi*S*rho*l).^2.*tmp;
    end
end