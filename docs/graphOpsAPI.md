# graphOps
### shortIntGraph
Convert the indices in a graph to 32-bit ints.  This takes less storage, but does not speed up much


```julia
shortIntGraph(a::SparseMatrixCSC{Tv,Ti<:Integer})
```

graphOps.jl:3



### floatGraph
Convert the nonzero entries in a graph to Float64


```julia
floatGraph(a::SparseMatrixCSC{Tv,Ti<:Integer})
```

graphOps.jl:6



### lap
Create a Laplacian matrix from an adjacency matrix. We might want to do this differently, say by enforcing symmetry


```julia
lap(a)
```

graphOps.jl:12



### unweight
Create a new graph in that is the same as the original, but with all edge weights 1


```julia
unweight{Tval,Tind}(ain::SparseMatrixCSC{Tval,Tind})
```

graphOps.jl:16



### unweight!
Change the weight of every edge in a to 1


```julia
unweight!{Tval,Tind}(a::SparseMatrixCSC{Tval,Tind})
```

graphOps.jl:26



### mapweight
Create a new graph that is the same as the original, but with f applied to each nonzero entry of a. For example, to make the weight of every edge uniform in [0,1], we could write

```julia
b = mapweight(a, x->rand(1)[1])
```


```julia
mapweight{Tval,Tind}(a::SparseMatrixCSC{Tval,Tind}, f)
```

graphOps.jl:40



### uniformWeight
Put a uniform [0,1] weight on every edge.  This is an example of how to use mapweight.


```julia
uniformWeight{Tval,Tind}(a::SparseMatrixCSC{Tval,Tind})
```

graphOps.jl:49



### uniformWeight!
Set the weight of every edge to 1


```julia
uniformWeight!(mat::SparseMatrixCSC{Tv,Ti<:Integer})
```

graphOps.jl:53



### productGraph
The Cartesian product of two graphs.  When applied to two paths, it gives a grid.


```julia
productGraph(a0::SparseMatrixCSC{Tv,Ti<:Integer}, a1::SparseMatrixCSC{Tv,Ti<:Integer})
```

graphOps.jl:58



### edgeVertexMat
The signed edge-vertex adjacency matrix


```julia
edgeVertexMat(mat::SparseMatrixCSC{Tv,Ti<:Integer})
```

graphOps.jl:67



### subsampleEdges
Create a new graph from the old, but keeping edge edge with probability `p`


```julia
subsampleEdges(a::SparseMatrixCSC{Float64,Int64}, p::Float64)
```

graphOps.jl:76



### twoLift
Creats a 2-lift of a.  `flip` is a boolean indicating which edges cross


```julia
twoLift(a)
twoLift(a, flip::AbstractArray{Bool,1})
twoLift(a, k::Integer)
```

graphOps.jl:99



### joinGraphs
Create a disjoint union of graphs a and b,  and then put k random edges between them


```julia
joinGraphs{Tval,Tind}(a::SparseMatrixCSC{Tval,Tind}, b::SparseMatrixCSC{Tval,Tind}, k::Integer)
```

graphOps.jl:111



### disjoin
Create a disjoint union of graphs a and b,   with no edges between them.


```julia
disjoin(a, b)
```

graphOps.jl:130



### plotGraph
Plots graph gr with coordinates (x,y)


```julia
plotGraph(gr, x, y)
plotGraph(gr, x, y, color)
```

graphOps.jl:136



### spectralDrawing
Computes spectral coordinates, and then uses plotGraph to draw


```julia
spectralDrawing(a)
```

graphOps.jl:173



### spectralCoords
Computes the spectral coordinates of a graph


```julia
spectralCoords(a)
```

graphOps.jl:181



### toUnitVector
Creates a unit vector of length n from a given set of integers, with weights based on the number of occurences


```julia
toUnitVector(a::Array{Int64,1}, n)
```

graphOps.jl:190



### diagmat
Returns the diagonal matrix(as a sparse matrix) of a graph


```julia
diagmat{Tv,Ti}(G::SparseMatrixCSC{Tv,Ti})
```

graphOps.jl:212



### generalizedNecklace
Constructs a generalized necklace graph starting with two graphs A and H. The resulting new graph will be constructed by expanding each vertex in H to an instance of A. k random edges will be generated between components. Thus, the resulting graph may have weighted edges.


```julia
generalizedNecklace{Tv,Ti}(A::SparseMatrixCSC{Tv,Ti}, H::SparseMatrixCSC{Tv,Ti<:Integer}, k::Int64)
```

graphOps.jl:232



