%
J = 0:20;
L = [0.1, 0.05, 0.025];

Y = zeros(length(L),length(J));
for i=1:length(L)
    l = L(i);
    Y(i,:) = exp(-(2*pi*l)^(-2))*besseli(J,(2*pi*l)^(-2));
end

figure
plot(J,Y,'o')

fid = fopen('pse_eigenvalues.dat','w');
colnames = 'j';
format = '%d';
for i=1:length(L)
    colnames = sprintf('%s\t%.3f',colnames,L(i));
    format = [format,'\t','%f'];
end
format = [format,'\n'];
fprintf(fid,'%s\n',colnames);
fprintf(fid,format,[J*2;Y]);
fclose(fid);