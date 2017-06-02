load fromJulia.mat

startTime = datestr(now);
st = Inf;
bt = Inf;
iter = Inf;
relres = Inf;
save toJulia st bt iter relres startTime

% run it once to warm it up ?
n0 = 1000;
a = rand(n0);
a = a + a';
a = a - diag(diag(a));
la0 = sparse(diag(sum(a)) - a);
b0 = randn(n0,1);
b0 = b0 - mean(b0);


inputType = 'laplacian';        % The input matrix A is a graph Laplacian
lamg    = Solvers.newSolver('lamg', 'randomSeed', 1);

setup   = lamg.setup('laplacian', la0);
[x, ~, ~, details] = lamg.solve(setup, b0, 'errorReductionTol', 0.1);


lamg    = Solvers.newSolver('lamg', 'randomSeed', 1);

%tic();
%setup   = lamg.setup('laplacian', la);
%bt = toc()
%lamg    = Solvers.newSolver('lamg', 'randomSeed', 1);

tic();
setup   = lamg.setup('laplacian', la);
bt = toc()


tic();

[x, ~, ~, details] = lamg.solve(setup, b, 'errorReductionTol', tol);

st = toc();

iter = length(details.stats.errorNormHistory);
relres = norm(la*x-b)/norm(b);

save toJulia st bt iter relres startTime

exit
