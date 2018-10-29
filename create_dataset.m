clear param

param.use_dbc = true; % ディリクレ境界条件を入れるか否か

param.sigma_eps = 0.1;

% % use se
% param.k.fh = @se_k;
% param.k.sigma = 1.0;
% param.k.l = 1.0;

% use matern
param.k.fh = @matern_k;
param.k.sigma = 1.0;
param.k.l = 1.0;
param.k.nu = 5;


param.domain=[0,10];
param.data_points = 50; % 訓練データ数
param.test_points = 1000; % テストデータ数
param.sampling_location = 'linear'; % 等間隔でサンプルを得る
% param.sampling_location = 'random'; % 一様乱数の位置でサンプルを得る

a = param.domain(1);
b = param.domain(2);

switch param.sampling_location
  case 'linear'
    x = linspace(a,b,param.data_points);
  case 'random'    
    x = a + rand(1,param.data_points)*(b-a);
end
x_test = linspace(a,b,param.test_points);
X = [x(:);x_test(:)];

if param.use_dbc
    bkup = param.sigma_eps;
    param.sigma_eps = 0; % 観測ノイズを0へ
    [mu,K] = std_gpr_fixed_hyps([a;b],[0;0],X,param);
    param.sigma_eps = bkup;
    F = mvnrnd(mu,K);
else
    X = [x(:);x_test(:)];
    K = param.k.fh(X,X,param);
    F = mvnrnd(zeros(1,n),K);
end

f = F(1:param.data_points);
y = f+param.sigma_eps*randn(size(f));
f_test = F(param.data_points+1:end);
y_test = f_test+param.sigma_eps*randn(size(f_test));

x=x(:); y=y(:); f=f(:);
x_test=x_test(:); f_test=f_test(:); y_test=y_test(:);
save('data.mat','x','f','y','x_test','f_test','y_test','param');

plot(x,y,'ro',x_test,f_test)