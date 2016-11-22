function IccAndCmg(n, maxtime, tol)
% function IccAndCmg(n, maxtime, tol)
%
% run the icc and cmg solvers on graphs of size n for a time up to
% maxtime, measured in seconds

%addpath ~/Laplacians/speedTests/
%addpath ~/Laplacians/matlab/
%init

if nargin < 3
    tol = 1e-6
end


opts.tol = tol;
opts.maxit = 100000;

t0 = now;

fn = ['IccAndCmg', num2str(n), '.csv'];

i = 0;

ns = [];
is = [];
nnzs = [];
iccBuild = [];
iccSolve = [];
cmgBuild = [];
cmgSolve = [];

while (now - t0 < maxtime/60/60/24)
    
    
    i = i + 1

    a = wtedChimera(n,i);
    
    ns = [ns;n];
    is = [is;i];
    nnzs = [nnzs;nnz(a)];
    
    la = lap(a);
    
    b = randn(n,1);
    b = b - mean(b);

    tic; f = iccSolver(la,[],opts); tb = toc;
    tic; x = f(b); ts = toc;
    
    iccBuild = [iccBuild ; tb];
    iccSolve = [iccSolve ; ts];
    
    tic; f = cmgSolver(la,[],opts); tb = toc;
    tic; x = f(b); ts = toc;
    
    cmgBuild = [cmgBuild ; tb];
    cmgSolve = [cmgSolve ; ts];
    
end

T = table(ns,is,nnzs,iccBuild,iccSolve,cmgBuild,cmgSolve)
writetable(T,fn)

