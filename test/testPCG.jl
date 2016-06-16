srand(1)
a = mapweight(grid2(20),x->1/(rand(1)[1]));
la = lap(a)
n = size(la)[1]
b = randn(n)
b = b - mean(b);

d = diag(la)
pre = (x)->(x ./ d)

x = pcg(la,b,pre,maxits=2000,maxtime=1,verbose=true);
x = Laplacians.pcgBLAS(la,b,pre,maxits=2000,maxtime=1,verbose=true);
x = Laplacians.pcgSlow(la,b,pre,maxits=2000,maxtime=1,verbose=true);
x = pcg(la,b,pre,maxits=2000,verbose=true);
x = Laplacians.pcgBLAS(la,b,pre,maxits=2000,verbose=true);
x = Laplacians.pcgSlow(la,b,pre,maxits=2000,verbose=true);
x = Laplacians.pcgBLAS(la,b,pre,tol=1e-2);
x = Laplacians.pcgSlow(la,b,pre,tol=1e-2);
x = pcg(la,b,pre,tol=1e-2);

x = cg(la,b,verbose=true,maxtime=1,maxits=1000);
x = Laplacians.cgBLAS(la,b,verbose=true,maxtime=1,maxits=1000);
x = Laplacians.cgSlow(la,b,verbose=true,maxtime=1,maxits=1000);
x = cg(la,b,verbose=true,maxits=2000);
x = Laplacians.cgBLAS(la,b,verbose=true,maxits=2000);
x = Laplacians.cgSlow(la,b,verbose=true,maxits=2000);
x = cg(la,b,verbose=true,tol=0.5);
x = Laplacians.cgBLAS(la,b,verbose=true,tol=0.5);
x = Laplacians.cgSlow(la,b,verbose=true,tol=0.5);




n = 100
a = wtedChimera(n,3)
la = lap(a)

f = cgSolver(la,verbose=false)
b = randn(n)
b = b - mean(b)
x = f(b)
norm(la*x-b)

x = f(b,maxits=15,verbose=true, maxtime=1)

d = diagm(1./diag(la))
pre = x -> d*x

x = pcg(la, b, diagm(diag(la)),maxits=10)
x = pcg(la, b, pre,maxits=10,verbose=true)
f = pcgSolver(la, diagm(diag(la)) ,maxits=10,verbose=true)
x = f(b)
f = pcgSolver(la,pre,maxits=10,verbose=false)
x = f(b,verbose=true, maxits=1000, maxtime = 2)
