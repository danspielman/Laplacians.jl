# flow
### maxflow
implementation of Dinic's algorithm. computes the maximum flow and min-cut in G between s and t. we consider the adjacency matrix to be the capacity matrix 


```julia
maxflow{Tv,Ti}(G::SparseMatrixCSC{Tv,Ti}, s::Int64, t::Int64)
```

flow.jl:4



