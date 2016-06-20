function la = lap(a)
% function la = lap(a)
%
% return the laplacian of adjacency matrix a
%

la = diag(sum(a)) - a;

