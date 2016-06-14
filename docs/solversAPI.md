# solvers
### lapWrapSolver
Takes a solver for solving nonsingular sdd systems, and returns a solver for solving Laplacian systems. The optional args tol and maxits are not necessarily taken by all solvers.  But, if they are, one can pass them here


```julia
lapWrapSolver(solver)
lapWrapSolver(solver, la::AbstractArray{T,N})
lapWrapSolver(solver, la::AbstractArray{T,N}, b)
```

solvers.jl:178



### augmentTree
Takes as input a tree and an adjacency matrix of a graph. It then computes the stretch of every edge of the graph wrt the tree.  It then adds back the k edges of highest stretch, and k edges sampled according to stretch


```julia
augmentTree{Tv,Ti}(tree::SparseMatrixCSC{Tv,Ti}, mat::SparseMatrixCSC{Tv,Ti}, k::Ti)
```

solvers.jl:208



### augTreePrecon
This is an augmented spanning tree preconditioner for diagonally dominant linear systems.  It takes as optional input a tree growing algorithm. The default is a randomized variant of Kruskal. It adds back 2sqrt(n) edges via augmentTree. With the right tree, it should never be too bad.


```julia
augTreePrecon{Tv,Ti}(ddmat::SparseMatrixCSC{Tv,Ti})
```

solvers.jl:258



### augTreeSolver
This is the solver that calls augTreePrecon


```julia
augTreeSolver{Tv,Ti}(ddmat::SparseMatrixCSC{Tv,Ti})
```

solvers.jl:280



### augTreeLapPrecon
A version of augTreePrecon specialized for Laplacians


```julia
augTreeLapPrecon{Tv,Ti}(ddmat::SparseMatrixCSC{Tv,Ti})
```

solvers.jl:291



### augTreeLapSolver
This is the solver that calls augTreeLapPrecon.  A solver specialized for Laplacians.


```julia
augTreeLapSolver{Tv,Ti}(ddmat::SparseMatrixCSC{Tv,Ti})
```

solvers.jl:313



