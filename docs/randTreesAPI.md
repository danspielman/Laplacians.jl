# randTrees
### randishKruskal
A heuristic for computing low-stretch spanning trees.  Where Kruskal's MST algorithm adds edges in order of weight, this algorithm adds them at random with probability proportional to their weight.


```julia
randishKruskal{Tv,Ti}(mat::SparseMatrixCSC{Tv,Ti})
```

randTrees.jl:11



### randishPrim
A heuristic for computing low-stretch spanning trees.  Where Prim's MST algorithm grows a cluster by always adding the edge on the boundary of maximum weight, this algorithm adds a boundary edge with probability proportional to its weight.


```julia
randishPrim{Tval,Tind}(mat::SparseMatrixCSC{Tval,Tind})
```

randTrees.jl:40



