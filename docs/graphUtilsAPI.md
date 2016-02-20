# graphUtils
### wdeg
Finds the weighted degree of a vertex in the graph 


```julia
wdeg{Tv,Ti}(mat::SparseMatrixCSC{Tv,Ti}, v::Ti)
```

graphUtils.jl:19



### setValue
Sets the value of a certain edge in a sparse graph; value can be 0 without the edges dissapearing 


```julia
setValue{Tv,Ti}(mat::SparseMatrixCSC{Tv,Ti}, v::Ti, i::Ti, a::Tv)
```

graphUtils.jl:29



### backIndices
Same as the above, but now the graph is in adjacency list form 

Computes the back indices in a graph in O(M+N). works if for every edge (u,v), (v,u) is also in the graph 


```julia
backIndices{Tv,Ti}(G::SparseMatrixCSC{Tv,Ti})
backIndices{Tv1,Tv2}(G::Array{Array{Tuple{Tv1,Tv2},1},1})
```

graphUtils.jl:35



### findEntries
Similar to findnz, but also returns 0 entries that have an edge in the sparse matrix 


```julia
findEntries{Tv,Ti}(G::SparseMatrixCSC{Tv,Ti})
```

graphUtils.jl:114



### compConductance
Returns the quality of the cut for a given graph and a given cut set s.   the result will be |outgoing edges| / min(|vertices in set|, |N - vertices in set|)


```julia
compConductance{Tv,Ti}(G::SparseMatrixCSC{Tv,Ti}, s::Array{Int64,1})
```

graphUtils.jl:136



### getVolume
Computes the volume of subset s in an unweighted graph G 


```julia
getVolume{Tv,Ti}(G::SparseMatrixCSC{Tv,Ti}, s::Array{Int64,1})
```

graphUtils.jl:149



### getObound
Computes the number of edges leaving s 


```julia
getObound{Tv,Ti}(G::SparseMatrixCSC{Tv,Ti}, s::Array{Int64,1})
```

graphUtils.jl:167



