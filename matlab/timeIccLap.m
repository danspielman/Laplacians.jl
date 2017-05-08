load fromJulia.mat

startTime = datestr(now);
st = Inf;
bt = Inf;
iter = Inf;
relres = Inf;
save toJulia st bt iter relres startTime

tic();
n = length(la);
lasub = la(1:(n-1),1:(n-1));
p = symrcm(lasub);
laperm = lasub(p,p);
L = ichol(laperm);

bt = toc();

tic();
bsub = b(p) - mean(b);
[xs,flag,relres,iter] = pcg(laperm, bsub, tol, maxits, L, L');
x = zeros(n,1);
x(p) = xs;
x = x - mean(x);
st = toc();

relres = norm(la*x-b)/norm(b)


save toJulia st bt iter relres startTime

exit
