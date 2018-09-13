%
S = -15:15;
L = [0.1, 0.05, 0.025];

Y = zeros(length(L),length(S));
for i=1:length(L)
    l = L(i);
    Y(i,:) = exp(-(2*pi*l)^(-2))*besseli(S,(2*pi*l)^(-2));
end

figure
stem(S,Y','o')

fid = fopen('pse_spectrum.dat','w');
colnames = 's';
format = '%d';
for i=1:length(L)
    colnames = sprintf('%s\t%.3f',colnames,L(i));
    format = [format,'\t','%f'];
end
format = [format,'\n'];
fprintf(fid,'%s\n',colnames);
fprintf(fid,format,[S;Y]);
fclose(fid);