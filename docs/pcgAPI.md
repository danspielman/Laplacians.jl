# pcg
### cg
`cg(mat, b; tol, maxits, maxtime, verbose)` solves a symmetric linear system. `tol` is set to 1e-6 by default, `maxits` defaults to Inf `maxtime` defaults to Inf.  It measures seconds. `verbose` defaults to false


```julia
cg(mat, b::Array{Float64,1})
cg(mat, b::Array{Float32,1})
cg(mat, b)
```

pcg.jl:45



### cgSolver
`cgSolver(mat; tol, maxits, maxtime, verbose)` creates a solver for a PSD system `mat`. The parameters are as described in cg.


```julia
cgSolver(mat)
```

pcg.jl:57



### pcg
`pcg(mat, b, pre; tol, maxits, maxtime, verbose)` solves a symmetric linear system using preconditioner `pre`. `pre` can be a function or a matrix.  If a matrix, a function to solve it is created with cholFact. `tol` is set to 1e-6 by default, `maxits` defaults to Inf `maxtime` defaults to Inf.  It measures seconds. `verbose` defaults to false


```julia
pcg(mat, b, pre::AbstractArray{T,N})
pcg(mat, b::Array{Float64,1}, pre::Function)
pcg(mat, b::Array{Float32,1}, pre::Function)
pcg(mat, b, pre::Function)
```

pcg.jl:63



### pcgSolver
`pcgSolver(mat, pre; tol, maxits, maxtime, verbose)` creates a solver for a PSD system using preconditioner `pre`. The parameters are as described in pcg.


```julia
pcgSolver(mat, pre)
```

pcg.jl:85



### pcgLapSolver
Create a solver that uses cg to solve Laplacian systems in mat. Specialized for the case when pre is a Laplacian matrix.  Fix the default parameters of the solver as given


```julia
pcgLapSolver(mat, pre)
```

pcg.jl:91



