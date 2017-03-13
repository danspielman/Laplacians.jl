[TOC]

# Solving linear equations in Laplacians and SDD matrices

The main purpose of this package is to experiment with the implementation of algorithms for solving systems of linear equations in Laplacian and symmetric, diagonally dominant, M-matrices (SDDM).

We recall that a matrix $ L $ is a _Laplacian_ matrix if:

* It is symmetric,
* its off-diagonal entries are non-positive, and
* all of its row sums are 0.

These conditions imply that the diagonal entries are non-negative, and that the matrix is singular.  So, we only hope to solve equations of the form  $ Lx = b $ when $b$ is in the span of the matrix.  When the graph of the nonzero entries of the matrix is connected, this is precisely when the sum of the entries in $ b $ is zero.  Laplacian matrices are always positive semidefinite.

A matrix $ M $ is a symmetric _M-matrix_ if:

* It is symmetric,
* its off-diagonal entries are non-positive, and
* it is positive definite.

A matrix symmetric $ M $ is _diagonally dominant_ if each of its diagonals is at least the sum of the absolute values of the off-diagonal entries in its row.  A Laplacians matrix is diagonally dominant.  A diagonally dominant matrix is always positive semidefinite.

A _SDDM_ matrix (symmetric, diagonally-dominant M-matrix) is a matrix that is both diagonally dominant and an M-matrix.  You may think of a SDDM matrix as a Laplacian plus a non-negative, non-zero, diagonal matrix.  However, this is only guaranteed to produce a SDDM matrix when the graph underlying the Laplacian is connected.

Laplacians.jl contains code for solving systems of linear equations in both Laplacian and SDDM matrices.  In fact, these problems are equivalent.  So, usually a solver for one type of system is implemented, and then wrapped to solve the other.
The same ideas can be used to solve systems of equations in SDD matrices (the off-diagonals can be positive or negative), but a wrapper for these has not yet been written.


## The Solver Interface

All of the SDDM solvers take the SDDM matrix as input.

_All of the Laplacian solvers take the adjacency matrix of the underlying graph as input._

To solve a system of linear equations, one first passes the matrix defining the system to a linear equation solving algorithm.  This will return a function that solves systems of linear equations in that matrix.  For example,

~~~julia
julia> n = 1000;
julia> a = wtedChimera(n);  # produces a graph, as a sparse adjacency matrix
julia> b = randn(n); 
julia> b = b - mean(b); # so there is a solution
julia> f = cholLap(a)
(::#79) (generic function with 1 method)
julia> x = f(b);
julia> la = lap(a);  # construct the Laplacian of a
julia> norm(la*x-b)
2.9565023548855584e-13
~~~

All of the solvers take the following keyword arguments.
This means that they are optional, and will be set to their default values if not specified.

* `tol` : the relative accuracy required: $ \| M x - b \| / \| b \| $.
* `maxits` : the maximum number of iterations, for iterative algorithms.
* `maxtime` : quit if it takes more than this many seconds.  Not all routines obey this, but they try.
* `verbose` : set to `true` to display diagnostic information.
* `pcgIts` : If the algorithm is iterative, this allows it to return the number of iterations it performed.  If `pcgIts` is an array of positive length, then its first entry is set to the number of iterations.  Where `verbose` prints this information, `pcgIts` allows it to be returned to other code.  To disable this set `pcgIts` to a zero length array, like `Int[]`.

Most of the solvers are iterative, exploiting the preconditioned conjugate gradient. These are the solvers for which `maxits`, `maxtime` and `pcgIts` make the most sense.  Some solvers, like Cholesky factorization, just ignore these parameters.

All of these parameters may be set in the call that constructs `f`.  They may then be over-ridden by again setting them in the call to `f`.
Let's see how this works when using the conjugate gradient.

~~~julia
julia> f = cgLapSolver(a, tol=1e-2, verbose=true)
(::f) (generic function with 1 method)
julia> x = f(b);
CG BLAS stopped after: 78 iterations with relative error 0.009590493139133275.
julia> norm(la*x-b)/norm(b)
0.00959049313913375

julia> pcgIts = [0]
1-element Array{Int64,1}:
 0
julia> x = f(b,verbose=false, pcgIts=pcgIts);
julia> pcgIts
1-element Array{Int64,1}:
 78

julia> x = f(b,verbose=true, maxits=50);
CG BLAS stopped after: 50 iterations with relative error 0.050483096216933886.

julia> x = f(b, tol=1e-4);
CG BLAS stopped after: 131 iterations with relative error 8.886882933346416e-5.
julia> norm(la*x-b)/norm(b)
8.886882933294668e-5
~~~  

 


For some experiments with solvers, including some of those below, look at the notebook Solvers.ipynb.

In the following, we document many of the solvers that have been implemented in this package.

## Cholesky Factorization

Cholesky factorization, the version of Gaussian Elimination for symmetric matrices, should be the first solver you try.  It will be very fast for matrices of dimension less than 1000, and for larger matrices coming from two-dimensional problems.

You can compute a cholesky factor directly with `cholfact`.  It does  more than just compute the factor, and it saves its result in a data structure that implements `\`.  It uses SuiteSparse by Tim Davis.

Here is an example of how you would use it to solve a general non-singular linear system.

~~~julia
a = grid2(5)
la = lap(a)
sddm = copy(la)
sddm[1,1] = sddm[1,1] + 1
F = cholfact(sddm)

n = size(a)[1]
b = randn(n)
x = F \ b
norm(sddm*x-b)

 	1.0598778281116327e-14
~~~

As `cholfact` does not satisfy our interface, we wrap it in a routine [`cholSDDM`](@ref) that does.

To solve systems in Laplacian matrices, use [`cholLap`](@ref).  Recall that this routine should be passed the adjacency matrix.


~~~julia
f = cholLap(a)
b = randn(n); 
b = b - mean(b);
norm(la*f(b) - b)
	2.0971536951312585e-15
~~~


## CG and PCG

We have written implementations of Conjugate Gradient (CG) and Preconditioned Conjugate Gradient (PCG) that satisfy the interface.
These routines use BLAS when possible, and slower routines when dealing with data types that BLAS cannot handle.  

~~~julia
srand(1)
n = 50
M = randn(n,n); M = M * M';
b = randn(n)
x = cg(M,b,maxits=100,verbose=true);
CG BLAS stopped after: 66 iterations with relative error 2.0166243927814765e-7.

bbig = convert(Array{BigFloat,1},b)
xbig = cg(M,bbig,maxits=100,tol=1e-30)
CG Slow stopped after: 50 iterations with relative error 2.18672511297479336887519117065525148757254642683072581090418060286711737398731e-38.

norm(M*xbig - bbig)
1.605742093628722039938504001423963138146137896744531914963345296279741402982296e-37
~~~

To create a function `f` that uses cg to solve systems in M, use [`cgSolver`](@ref).  For Laplacians, use [`cgLapSolver`](@ref).

~~~julia
julia> n = 1000;
julia> a = wtedChimera(n,1);
julia> f = Laplacians.cgLapSolver(a,maxits=100);

julia> b = randn(n);
julia> b = b - mean(b);
julia> x = f(b,verbose=true);
CG BLAS stopped after: 100 iterations with relative error 0.012102058751548373.


julia> la = lap(a);
julia> sddm = copy(la);
julia> sddm = sddm + spdiagm(rand(n)/100);
julia> g = cgSolver(sddm,verbose=true)
(::f) (generic function with 1 method)

julia> x = g(b);
CG BLAS stopped after: 253 iterations with relative error 7.860172210007891e-7.
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
prec(x) = x ./ d
@time x = pcg(la,b,prec,maxtime=1,tol=1e-2,verbose=true);

PCG BLAS stopped at maxtime.
PCG BLAS stopped after: 530 iterations with relative error 0.07732478003311881.
  1.007756 seconds (10.32 k allocations: 648.525 MB, 9.69% gc time)
 
@time x = pcg(la,b,prec,maxtime=3,tol=1e-2,verbose=true);
PCG BLAS stopped after: 1019 iterations with relative error 0.009984013184429813.
  2.086828 seconds (19.57 k allocations: 1.216 GB, 9.92% gc time) 
~~~

Without the preconditioner, CG takes much longer on this example.

~~~julia
@time x = cg(la,b,tol=1e-2,maxtime=10,verbose=true);

CG BLAS stopped at maxtime.
CG BLAS stopped after: 8879 iterations with relative error 0.054355534834831624.
 10.001998 seconds (97.91 k allocations: 2.649 GB, 4.48% gc time)
~~~

[`pcgSolver`](@ref) creates a function that uses the preconditioner to solve systems in the matrix.

~~~julia
f = pcgSolver(la,prec)
@time x = f(b,maxtime=3,tol=1e-2,verbose=true);
PCG BLAS stopped after: 1019 iterations with relative error 0.009984013184429813.
  1.892217 seconds (19.58 k allocations: 1.216 GB, 9.47% gc time)
~~~

[`pcgLapSolver`](@ref) uses the Laplacian of one matrix as a preconditioner for the first.  It solves systems of linear equations in the preconditioner by Cholesky factorization.  It performs the Cholesky factorization when [`pcgLapSolver`](@ref) is called.  This is why we do the work of creating `f` only once.  Here is an example using a Low-Stretch Spanning Tree preconditioner.

~~~julia
@time t = akpw(a)
  0.210467 seconds (1.43 M allocations: 91.226 MB, 19.23% gc time)

@time f = pcgLapSolver(a,t)
  0.160210 seconds (288 allocations: 28.076 MB, 72.28% gc time)

@time x = f(b,maxtime=3,tol=1e-2,verbose=true);
PCG BLAS stopped after: 260 iterations with relative error 0.009864463201800925.
  1.014897 seconds (28.02 k allocations: 1.008 GB, 9.81% gc time)
~~~



## Low-Stretch Spanning Trees

In order to make preconditioners, we will want low-stretch spanning trees.  We do not yet have any code that is guaranteed to produce these.  Instead, we supply three heuristics: [`akpw`](@ref) which is inspired by the algorith of Alon, Karp, Peleg and West, and  randomized versions of Prim and Kruskal's algorithm.
[`randishKruskal`](@ref) samples the remaining edges with probability proportional to their weight.  [`randishPrim`](@ref) samples edges on the boundary while using the same rule.  We recommend using [`akpw`](@ref).

See [Low Stretch Spanning Trees](LSST.md) to learn more about these.

## Augmented Spanning Tree Preconditioners

These are obtained by constructing a spanning tree of a graph, and then adding back some more edges from the graph.  The tree should have low stretch.  The edges to add back are chosen at random with probabilities proportional to their stretches.

These are implemented in the routines

* [`augTreeSolver`](@ref), for SDDM matrices
* [`augTreeLapSolver`](@ref)
* [`augTreePrecon`](@ref)
* [`augmentTree`](@ref)


## The solvers of Koutis, Miller and Peng.

Solvers inspired by the algorithm from "Approaching optimality for solving SDD systems" by Koutis, Miller, and Peng, <i>SIAM Journal on Computing</i>, 2014.

* [`KMPLapSolver`](@ref)
* [`KMPSDDMSolver`](@ref)

## Sampling Solvers of Kyng and Sachdeva

These are inspired by the paper "Approximate Gaussian Elimination for Laplacians: Fast, Sparse, and Simple" by Rasmus Kyng and Sushant Sachdeva, FOCS 2016. 

* [`samplingSDDMSolver`](@ref)
* [`samplingLapSolver`](@ref)

## Algebraic Multigrid

This is an interface to the algebraic multigrid solvers from the PyAMG package.

* [`AMGLapSolver`](@ref)
* [`AMGSolver`](@ref), for SDDM systems.

## Solvers from Matlab

The [MATLAB.jl](https://github.com/JuliaInterop/MATLAB.jl) package allows Julia to call routines from Matlab, provided you have Matlab installed.  It does this in a very efficient fashion: it starts up the Matlab process when you type `using MATLAB`, and then communicates with it.  So, we have wrapped some solvers from Matlab so that they obey the same interface.

These are not part of the Laplacians module, but are included in the package under `src/matlabSolvers.jl`.  To include them, type

~~~julia
include(string(Pkg.dir("Laplacians") , "/src/matlabSolvers.jl"))
~~~

We provide the docstrings for these here.

### Incomplete Cholesky Factorizations

These use the no-fill incomplete Cholesky factorizations implemented in Matlab.  They first order the vertices by the `symrcm` ordering.


The solvers are:

* `f = matlab_ichol_sddm(sddm; tol, maxtime, maxits, pctIts, verbose)`
* `f = matlab_ichol_lap(A; tol, maxtime, maxits, pctIts, verbose)`

A routine that just wraps the function that solves equations in the preconditioner is provided as well:

* `f = matlab_ichol(sddm)`

### Koutis's Combinatorial Multigrid (CMG)

You must have installed Yiannis Koutis's [Combinatorial Multigrid Code](http://www.cs.cmu.edu/~jkoutis/cmg.html), and it must be on Matlab's default path.  As this code returns a function rather than a preconditioner, it would be inefficient to make it use our PCG code and satisfy our interface.  So, it does not.

* `x = matlabCmgSolver(mat, b; tol::Real=1e-6, maxits=10000)`

The matrix `mat` can either be SDDM or a Laplacian.  This solves the system in `b`.
 