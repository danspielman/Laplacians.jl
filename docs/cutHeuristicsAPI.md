# cutHeuristics
### refineCut
```
Modify a cluster by adding or removing vertices by picking at each step 
the vertex that has the maximum value of (Deg_external - Deg_Internal).
Each vertex can be added in/removed only once.
```


```julia
refineCut{Tv,Ti}(G::SparseMatrixCSC{Tv,Ti}, s::Array{Int64,1})
```

cutHeuristics.jl:9



### dumb
```
Modify a cluster by passing through all the vertices exactly once and 
adding/removing them based on the value of (Deg_external - Deg_Internal).
```


```julia
dumb{Tv,Ti}(G::SparseMatrixCSC{Tv,Ti}, s::Array{Int64,1})
```

cutHeuristics.jl:106



