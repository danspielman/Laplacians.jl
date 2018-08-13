Random.seed!(1)
a = mapweight(grid2(20),x->1/(rand(1)[1]));
la = lap(a)
n = size(la)[1]
b = randn(n)
b = b .- mean(b);

d = Array(diag(la))
pre = (x)->(x ./ d)

x = pcg(la,b,pre,maxits=2000,maxtime=1,verbose=true);
x = Laplacians.pcgBLAS(la,b,pre,maxits=2000,maxtime=1,verbose=true);
x = Laplacians.pcgSlow(la,b,pre,maxits=2000,maxtime=1,verbose=true);
x = pcg(la,b,pre,maxits=2000,verbose=true);
x = Laplacians.pcgBLAS(la,b,pre,maxits=2000,verbose=true);
x = Laplacians.pcgBLAS(la,b,pre,maxits=2000,verbose=true,maxtime=1e-10);
x = Laplacians.pcgSlow(la,b,pre,maxits=2000,verbose=true);
x = Laplacians.pcgSlow(la,b,pre,maxits=2000,verbose=true,maxtime=1e-10);
x = Laplacians.pcgBLAS(la,b,pre,tol=1e-2);
x = Laplacians.pcgSlow(la,b,pre,tol=1e-2);
x = pcg(la,b,pre,tol=1e-2);

x = cg(la,b,verbose=true,maxtime=1,maxits=1000);
x = Laplacians.cgBLAS(la,b,verbose=true,maxtime=1,maxits=1000);
x = Laplacians.cgBLAS(la,b,verbose=true,maxtime=1e-10,maxits=1000);
x = Laplacians.cgSlow(la,b,verbose=true,maxtime=1,maxits=1000);
x = Laplacians.cgSlow(la,b,verbose=true,maxtime=1e-10,maxits=1000);
x = cg(la,b,verbose=true,maxits=2000);
x = Laplacians.cgBLAS(la,b,verbose=true,maxits=2000);
x = Laplacians.cgSlow(la,b,verbose=true,maxits=2000);
x = cg(la,b,verbose=true,tol=0.5);
x = Laplacians.cgBLAS(la,b,verbose=true,tol=0.5);
x = Laplacians.cgSlow(la,b,verbose=true,tol=0.5);

bbig = convert(Array{BigFloat,1},b);
f = cgSolver(la)
x = f(bbig);

x = pcg(la,bbig,pre,maxits=200,maxtime=1,verbose=true);
x = cg(la,bbig,maxits=200,maxtime=1,verbose=true);


n = 100
a = wted_chimera(n,3)
la = lap(a)

f = cgSolver(la,verbose=false)
b = randn(n)
b = b .- mean(b)
x = f(b)
norm(la*x-b)

f(zeros(100))

x = f(b,maxits=15,verbose=true, maxtime=1)

d = Diagonal(1 ./ diag(la))
pre = x -> Array(d*x)

x = pcg(la, b, sparse(Diagonal(diag(la))),maxits=10)
x = pcg(la, b, pre, maxits=10, verbose=true)
f = pcgSolver(la, spdiagm(0=>diag(la)) ,maxits=10,verbose=true)
x = f(b)
f = pcgSolver(la,pre,maxits=10,verbose=false)
x = f(b,verbose=true, maxits=1000, maxtime = 2)

f(zeros(100))

t = akpw(a)
lt = lap(t)

f = pcgLapSolver(a, t)
pcgIts = [1]
x = f(b,verbose=true, pcgIts=pcgIts)
bbig = convert(Array{BigFloat,1},b);
x = f(bbig,verbose=true, pcgIts=pcgIts)

sdda = copy(la)
sdda[1,1] += 1
sddt = copy(lt)
sddt[1,1] += 1

f = pcgSolver(sdda, sddt)
x = f(b,verbose=true, pcgIts=pcgIts)
x = f(bbig,verbose=true, pcgIts=pcgIts)


