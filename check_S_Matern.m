% Maternカーネルのスペクトル密度とt分布

m = 10;
x = linspace(0,.5*m,51);

Nu = [1,10];
L = [.05,.1];

figure; hold on
for i=1:length(L)
    l = L(i);
    for j=1:length(Nu)
        nu = Nu(j);
        
        S = S_Matern(x,nu,l);
        tau = 2*nu;
        sigma = 1/(2*pi*l);
        td = makedist('tLocationScale','mu',0,'sigma',sigma,'nu',tau);
        J = .5*(1:m); % 周期は 1/2,1,3/2,...
        p = pdf(td,J);
        plot(x,S,J,p,'o')
    end
end

% function S = S_Matern(j,nu,l)
%     S = 2^(nu+1)*nu^nu / (beta(nu,.5)*l^(2*nu)) ...
%         * (2*nu/l^2+4*pi^2*j.^2).^-(nu+.5);
% end
function S = S_Matern(j,nu,l)
    S = 2*sqrt(pi)*gamma(nu+.5)*(2*nu)^nu / (gamma(nu)*l^(2*nu)) ...
        * (2*nu/l^2+4*pi^2*j.^2).^(-nu-.5);
end