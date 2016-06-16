srand(1)
a = mapweight(grid2(200),x->1/(rand(1)[1]));
la = lap(a)
n = size(la)[1]
b = randn(n)
b = b - mean(b);

d = diag(la)
pre(x) = x ./ d

x = pcg(la,b,pre,maxits=2000,maxtime=1,verbose=true);
x = pcgBLAS(la,b,pre,maxits=2000,maxtime=1,verbose=true);
x = pcgSlow(la,b,pre,maxits=2000,maxtime=1,verbose=true);
x = pcg(la,b,pre,maxits=2000,verbose=true);
x = pcgBLAS(la,b,pre,maxits=2000,verbose=true);
x = pcgSlow(la,b,pre,maxits=2000,verbose=true);
x = pcgBLAS(la,b,pre,tol=1e-2);
x = pcgSlow(la,b,pre,tol=1e-2);
x = pcg(la,b,pre,tol=1e-2);

x = cg(la,b,verbose=true,maxtime=1,maxits=1000);
x = cgBLAS(la,b,verbose=true,maxtime=1,maxits=1000);
x = cgSlow(la,b,verbose=true,maxtime=1,maxits=1000);
x = cg(la,b,verbose=true,maxits=2000);
x = cgBLAS(la,b,verbose=true,maxits=2000);
x = cgSlow(la,b,verbose=true,maxits=2000);
x = cg(la,b,verbose=true,tol=0.5);
x = cgBLAS(la,b,verbose=true,tol=0.5);
x = cgSlow(la,b,verbose=true,tol=0.5);




