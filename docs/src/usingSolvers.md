
# Solving linear equations in Laplacians


For some experiments with solvers, including some of those below, look at the notebook Solvers.ipynb.

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

Our implementation requires 2 arguments: the matrix and the right-hand vector.  It's optional arguments are the tolerance `tol`, the maximum number of iterations, `maxits`, and the maximum number of seconds, `maxtime`.  It has been written to use BLAS when possible, and slower routines when dealing with data types that BLAS cannot handle.  For extra information, set `verbose` to true.  Here are examples.

~~~julia
srand(1)
n = 50
a = randn(n,n); a = a * a';
b = randn(n)
x = cg(a,b,maxits=100,verbose=true);
	CG stopped after: 66 iterations with relative error 4.901070762774203e-7.
norm(a*x - b)/norm(b)
	4.901071329720876e-7

bbig = convert(Array{BigFloat,1},b)
xbig = cg(a,bbig,maxits=100)
	CG stopped after: 50 iterations with relative error 2.18672511297479336887519117065525148757254642683072581090418060286711737398731e-38.
norm(a*xbig - bbig)
	2.186725112974793368875191170655251487433448469718749360094794057951719231569009e-38
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
srand(1)
a = mapweight(grid2(200),x->1/(rand(1)[1]));
la = lap(a)
n = size(la)[1]
b = randn(n)
b = b - mean(b);

d = diag(la)
pre(x) = x ./ d
@time x = pcg(la,b,pre,maxits=2000,verbose=true);

	PCG stopped after: 2000 iterations with relative error 1.535193885968155e-5.
	  2.719655 seconds (42.41 k allocations: 1.194 GB, 6.21% gc time)
~~~

If our target is just low error, and we are willing to allow many iterations, here's how cg and pcg compare on this example.

~~~julia
@time x = pcg(la,b,pre,tol=1e-1,maxits=10^5,verbose=true);
	PCG stopped after: 452 iterations with relative error 0.09944363573171432.
 	0.649586 seconds (9.67 k allocations: 277.034 MB, 5.74% gc time) 
norm(la*x - b)/norm(b)
	0.09944363573090279

@time x = cg(la,b,tol=1e-1,maxits=10^5);
	4.526119 seconds (16.03 k allocations: 1.195 GB, 4.16% gc time)
  
norm(la*x - b)/norm(b)
	0.0981610595725004  
~~~

## Low-Stretch Spanning Trees

In order to make preconditioners, we will want low-stretch spanning trees.  We do not yet have any code that is guaranteed to produce these.  Instead, we supply three heuristics: `akpw` which is inspired by the algorith of Alon, Karp, Peleg and West, and  randomized versions of Prim and Kruskall's algorithm.
`randishKruskall` samples the remaining edges with probability proportional to their weight.  `randishPrim` samples edges on the boundary while using the same rule.

See [Low Stretch Spanning Trees](LSST.md)

Let's try using a low-stretch spanning tree as a preconditioner for that last example.  First, we compute a low stretch tree and, for fun, check its average stretch.

~~~julia
@time tree = akpw(a);
	0.195370 seconds (1.07 M allocations: 86.283 MB, 5.84% gc time)
aveStretch = sum(compStretches(tree,a))/nnz(a)
	6.527424129533183
~~~

To use it as a preconditioner, we must construct a function that inverts linear systems in its Laplacian:

~~~julia
ltree = lap(tree)
@time ltreeSolver = lapChol(ltree);
	0.045963 seconds (92 allocations: 17.703 MB, 9.51% gc time)
norm(ltree*ltreeSolver(b) - b)
	8.209228285591316e-8
~~~

Now, let's see how well it works as a preconditioner.

~~~julia
@time x = pcg(la,b,ltreeSolver,tol=1e-1,maxits=10^5,verbose=true);
	PCG stopped after: 153 iterations with relative error 0.09900131423570005.
	0.658006 seconds (12.27 k allocations: 654.998 MB, 11.53% gc time)
~~~

That's a lot like using the diagonal preconditioner.  However, we see that the tree performs fewer iterations, and is faster when we seek higher accuracy.

~~~julia
@time x = pcg(la,b,pre,tol=1e-6,maxtime=10,verbose=true);
	PCG stopped after: 2435 iterations with relative error 9.92454258869586e-7.
 	3.431920 seconds (51.35 k allocations: 1.454 GB, 6.16% gc time)

@time x = pcg(la,b,ltreeSolver,tol=1e-6,maxtime=10,verbose=true);
	PCG stopped after: 605 iterations with relative error 9.78414046170984e-7.
 	2.510224 seconds (48.01 k allocations: 2.527 GB, 11.68% gc time)
~~~


## Various Linear System Solvers

We have a blend of solvers geared towards solving symmetric linear systems. Some of these solvers, that have
historically showcased good performance, have been grouped into two sections: SDDSolvers and LapSolvers.

~~~julia
SDDSolvers
5-element Array{Function,1}:
  Laplacians.augTreeSolver    
  Laplacians.KMPSDDSolver     
  Laplacians.hybridSDDSolver  
  Laplacians.samplingSDDSolver
  Laplacians.AMGSolver

LapSolvers
5-element Array{Function,1}:
  Laplacians.augTreeLapSolver 
  Laplacians.KMPLapSolver     
  Laplacians.hybridLapSolver  
  Laplacians.samplingLapSolver
  Laplacians.AMGLapSolver 
~~~

The SDD solvers take in a SDD matrix. They can also take in the usual tol, maxits and maxtime parameters. The
Laplacian solvers, in contrast, are fed in an adjacency matrix. They are also subject to tol, maxits and maxtime.
Following are examples of how to run these solvers on a given linear system.

~~~julia
a = mapweight(grid2(100),x->1/(rand(1)[1]));
n = a.n
la = lap(a); 
dval = zeros(n); dval[1] = dval[n] = 1e-3;
sdd = la + spdiagm(dval);
b = randn(n); b = b - mean(b);
~~~

Now, to run the SDD solvers.

~~~julia
for solver in SDDSolvers
    println("Solver ", solver)
    @time f = solver(sdd, maxits=1000, maxtime=10, tol=1e-4, verbose=false)
    @time x = f(b);
    println("Relative norm: ", norm(sdd * x - b) / norm(b), "\n")
end

Solver Laplacians.augTreeSolver
  0.067958 seconds (235.22 k allocations: 35.807 MB, 8.41% gc time)
  0.060869 seconds (2.72 k allocations: 70.339 MB, 17.35% gc time)
Relative norm: 9.250884271915595e-5

Solver Laplacians.KMPSDDSolver
  0.063388 seconds (254.46 k allocations: 44.168 MB, 14.65% gc time)
  0.102330 seconds (8.82 k allocations: 133.219 MB, 18.42% gc time)
Relative norm: 9.355785387633888e-5

Solver Laplacians.hybridSDDSolver
  0.064385 seconds (264.54 k allocations: 41.390 MB, 12.33% gc time)
  0.379060 seconds (2.07 M allocations: 266.961 MB, 9.90% gc time)
Relative norm: 9.832304243953024e-5

Solver Laplacians.samplingSDDSolver
  0.135522 seconds (706.74 k allocations: 101.318 MB, 13.69% gc time)
  0.179568 seconds (3.14 M allocations: 144.395 MB, 13.88% gc time)
Relative norm: 9.95716064856691e-5

Solver Laplacians.AMGSolver
  0.030372 seconds (166 allocations: 860.125 KB)
  0.116332 seconds (1.55 k allocations: 4.400 MB)
Relative norm: 9.519986283623698e-5
~~~

Next, the Laplacian solvers.

~~~julia
for solver in LapSolvers
    println("Solver ", solver)
    @time f = solver(a, maxits=1000, maxtime=10, tol=1e-4, verbose=false)
    @time x = f(b);
    println("Relative norm: ", norm(sdd * x - b) / norm(b), "\n")
end

Solver Laplacians.augTreeLapSolver
  0.072029 seconds (235.32 k allocations: 38.405 MB, 12.87% gc time)
  0.088307 seconds (6.51 k allocations: 98.563 MB, 18.56% gc time)
Relative norm: 0.00015469840172506333

Solver Laplacians.KMPLapSolver
  0.054602 seconds (254.18 k allocations: 39.445 MB, 15.07% gc time)
  0.107957 seconds (8.38 k allocations: 126.549 MB, 21.16% gc time)
Relative norm: 0.0001558714652838137

Solver Laplacians.hybridLapSolver
  0.053575 seconds (267.23 k allocations: 36.942 MB, 6.97% gc time)
  0.407590 seconds (3.02 M allocations: 286.329 MB, 10.58% gc time)
Relative norm: 0.00015420023270231535

Solver Laplacians.samplingLapSolver
  0.119625 seconds (706.35 k allocations: 90.641 MB, 21.37% gc time)
  0.218231 seconds (3.52 M allocations: 161.625 MB, 14.20% gc time)
Relative norm: 0.00015470313398451233

Solver Laplacians.AMGLapSolver
  0.043749 seconds (235 allocations: 2.744 MB)
  0.179676 seconds (1.57 k allocations: 4.402 MB)
Relative norm: 0.00014173955766018617
~~~
