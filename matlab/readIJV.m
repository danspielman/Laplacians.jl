function a = readIJV(filename)
% function a = readIJV(filename)
%
% read a comma-separated IJV from file,
% and return a sparse matrix

dl = dlmread(filename);
a = sparse(dl(:,1),dl(:,2),dl(:,3));
n = max(size(a));
a(n,n) = 0;
a = a + a';
