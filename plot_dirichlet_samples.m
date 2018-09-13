a=0; b=1;
param.domain=[a,b];
param.k.l = 0.1;
param.k.sigma = 1;

x = linspace(a,b,101);
C = se_k(x,x,param);
W = win_cos(x,x,param);
C = W.*C;

F = mvnrnd(zeros(size(x)),C,5);
plot(x,F')

fid = fopen('samples_coswin_x_se.dat','w');
fprintf(fid,'x\tf1\tf2\tf3\tf4\tf5\n');
fprintf(fid,'%f\t%f\t%f\t%f\t%f\t%f\n',[x;F]);
fclose(fid);
