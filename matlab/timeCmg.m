load fromJulia.mat

startTime = datestr(now);
st = Inf;
bt = Inf;
iter = Inf;
relres = Inf;
save toJulia st bt iter relres startTime

tic();
pfun = cmg_sdd(la);
bt = toc();

tic();
[x,flag,relres,iter] = pcg(la, b, tol, maxits, pfun);
st = toc();

save toJulia st bt iter relres startTime

exit
