# pcg
### cg
`cg(mat, b; tol, maxits, maxtime, verbose)` solves a symmetric linear system. `tol` is set to 1e-6 by default, `maxits` defaults to Inf `maxtime` defaults to Inf.  It measures seconds. `verbose` defaults to false


```julia
cg(mat, b::Array{Float64,1})
cg(mat, b::Array{Float32,1})
cg(mat, b)
```

pcg.jl:31



### pcg
`pcg(mat, b, pre; tol, maxits, maxtime, verbose)` solves a symmetric linear system using preconditioner `pre`. `pre` should be a function `tol` is set to 1e-6 by default, `maxits` defaults to Inf `maxtime` defaults to Inf.  It measures seconds. `verbose` defaults to false


```julia
pcg(mat, b::Array{Float64,1}, pre)
pcg(mat, b::Array{Float32,1}, pre)
pcg(mat, b, pre)
```

pcg.jl:44



