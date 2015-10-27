[TOC]

# Solving linear equations in Laplacians

Right now, our solver code is in `solvers.jl`, but not included in yinsGraph.  So, you should include this directly.  Implementations of cg and pcg have been automatically included in yinsGraph.  They are in the file `pcg.jl`

For some experiments with solvers, including some of those below, look at the notebook [Solvers.ipynb](Solvers.ipynb).

## Direct Solvers

You can compute a cholesky factor directly with `cholfact`.  It does  more than just compute the factor, and it saves its result in a data structure that implements `\`.  It uses SuiteSparse by Tim Davis.

Here is an example of how you would use it to solve a general non-singular linear system.

~~~julia
a = grid2(5)
la = lap(a)
la[1,1] = la[1,1] + 1
F = cholfact(la)

n = size(a)[1]
b = randn(n)
x = F \ b
norm(la*x-b)

 	1.0598778281116327e-14
~~~

Laplacians, however, are singular.  So, we need to wrap the solver inside a routine that compensates for this.

~~~julia
la = lap(a)
f = lapWrapSolver(cholfact,la)
b = randn(n); b = b - mean(b);
norm(la*f(b) - b)
	2.0971536951312585e-15
~~~

Here are two other ways of using the wrapper:

~~~julia
lapChol = lapWrapSolver(cholfact)
f = lapChol(la)
b = randn(n); 
b = b - mean(b);
norm(la*f(b) - b)
	2.6924696662484416e-15
	
x = lapWrapSolver(cholfact,la,b)
norm(la*x - b)
	2.6924696662484416e-15
~~~

## Iterative Solvers
The first, of course, is the Conjugate Gradient (cg).

Our implementation requires 2 arguments: the matrix and the right-hand vector.  It's optional arguments are the tolerance `tol` and the maximum number of iterations, `maxits`.  It has been written to use BLAS when possible, and slower routines when dealing with data types that BLAS cannot handle.  Here are examples.

~~~julia
n = 50
a = randn(n,n); a = a * a';
b = randn(n)
x = cg(a,b,maxits=100)
norm(a*x - b)
	1.2191649497921835e-6
	
bbig = convert(Array{BigFloat,1},b)
xbig = cg(a,bbig,maxits=100)
norm(a*xbig - bbig)
	1.494919244242202629856363570306545126541716514824419323325986374186529786019681e-33
~~~

As a sanity check, we do two speed tests against Matlab.

~~~julia
la = lap(grid2(200))
n = size(la)[1]
b = randn(n)
b = b - mean(b);
@time x = cg(la,b,maxits=1000)
	0.813791 seconds (2.77 k allocations: 211.550 MB, 3.56% gc time)

norm(la*x-b)
	0.0001900620047823064
~~~

And, in Matlab:

~~~matlab
>> a = grid2(200);
>> la = lap(a);
>> b = randn(length(a),1); b = b - mean(b);
>> tic; x = pcg(la,b,[],1000); toc
pcg converged at iteration 688 to a solution with relative residual 9.8e-07.
Elapsed time is 1.244917 seconds.
>> norm(la*x-b)

ans =

   1.9730e-04
~~~

PCG also takes as input a preconditioner.  This should be a function.  Here is an example of how one might construct and use a diagonal preonditioner.  To motivate this, I will use a grid with highly varying weights on edges.

~~~julia
a = mapweight(grid2(200),x->1/(rand(1)[1]));
la = lap(a)
n = size(la)[1]
b = randn(n)
b = b - mean(b);

d = diag(la)
pre(x) = x ./ d
@time x = pcg(la,b,pre,maxits=2000)
	3.322035 seconds (42.21 k allocations: 1.194 GB, 5.11% gc time)
norm(la*x-b)
	0.008508746034886803
~~~

If our target is just low error, and we are willing to allow many iterations, here's how cg and pcg compare on this example.

~~~julia
@time x = pcg(la,b,pre,tol=1e-1,maxits=10^5)
	0.747042 seconds (9.65 k allocations: 275.819 MB, 4.87% gc time)
norm(la*x-b)
	19.840756251253442
	
@time x = cg(la,b,tol=1e-1,maxits=10^5)
	6.509665 seconds (22.55 k allocations: 1.680 GB, 3.68% gc time)
norm(la*x-b)
	19.222483530605043
~~~

## Low-Stretch Spanning Trees

In order to make preconditioners, we will want low-stretch spanning trees.  We do not yet have any code in Julia that is guaranteed to produce these.  Instead, for now, we have two routines that can be thought of as randomized versions of Prim and Kruskall's algorithm.
`randishKruskall` samples the remaining edges with probability proportional to their weight.  `randishPrim` samples edges on the boundary while using the same rule.

Both use a data structure called `Sampler` that allows you to store integers with real values, and to sample according to those real values.

We also have code for computing the stretches.
Here are some examples.

~~~julia
a = grid2(1000)
t = randishKruskal(a);
st = compStretches(t,a);
sum(st)/nnz(a)
	43.410262262262265
	
t = randishPrim(a);
st = compStretches(t,a);
sum(st)/nnz(a)
	33.14477077077077

~~~

## Augmented spanning tree preconditioners

Here is code that will invoke one.  
It is designed for positive definite systems.  So, let's give it one.
Right now, it is using a randomized version of a MST.  There is no real reason to think that this should work.


~~~julia
a = mapweight(grid2(1000),x->1/(rand(1)[1]));
la = lap(a)
n = size(la)[1]
la[1,1] = la[1,1] + 1
@time F = augTreeSolver(la,tol=1e-1,maxits=1000)
	6.529052 seconds (4.00 M allocations: 1.858 GB, 15.34% gc time)

b = randn(n)
@time x = F(b)
	29.058915 seconds (9.74 k allocations: 23.209 GB, 6.84% gc time)

norm(la*x - b)
	99.74452367765869

# Now, let's contrast with using CG

@time y = cg(la,b,tol=1e-1,maxits=1000)
	28.719631 seconds (4.01 k allocations: 7.473 GB, 3.74% gc time)

norm(la*y-b)
	3243.6014713600766

~~~
That was not too impressive.  We will have to investigate.  By default, it presently uses randishKruskal.  Let's try randishPrim.  You can pass the treeAlg as a parameter.

~~~julia
@time F = augTreeSolver(la,tol=1e-1,maxits=1000,treeAlg=randishPrim);
	6.319489 seconds (4.00 M allocations: 2.030 GB, 18.81% gc time)
	
b = randn(n)
@time x = F(b)
	29.503484 seconds (9.76 k allocations: 23.268 GB, 7.31% gc time)

norm(la*x - b)	
	99.29610874176991
~~~


To solve systems in a Laplacian, we could wrap it.

~~~julia
n = 40000
la = lap(randRegular(n,3))
f = lapWrapSolver(augTreeSolver,la,tol=1e-6,maxits=1000)
b = randn(n); b = b - mean(b)
x = f(b)
norm(la*x-b)
	0.00019304778073388
~~~

As you can see, lapWrapSolver can pass tol and maxits arguments to its solver, if they are given to it.




