# localClustering
### localImprove
localImprove{Tv,Ti}(G::SparseMatrixCSC{Tv,Ti}, A::Array{Int64,1}; epsSigma=-1.0, err=1e-10, maxSize = max(G.n, G.m)

The LocalImprove function, from the Orrechia-Zhu paper 

G is the given graph, A is the initial set    epsSigma is a measure of the quality of the returning set (the smaller the better)   err is the numerical error considered throughout the algorithm   maxSize is the maximum allowed size for the flow graph at any iteration of the algorithm


```julia
localImprove{Tv,Ti}(G::SparseMatrixCSC{Tv,Ti}, A::Array{Int64,1})
```

localClustering.jl:16



### prn
prn{Tv, Ti}(G::SparseMatrixCSC{Tv,Ti}, s::Array{Int64,1}, phi::Float64, b::Int64)

The PageRank-Nibble cutting algorithm from the Anderson/Chung/Lang paper

s is a set of starting vertices, phi is a constant in (0, 1], and b is an integer in [1, [log m]]

phi is a bound on the quality of the conductance of the cut - the smaller the phi, the higher the quality   b is used to handle precision throughout the algorithm - the higher the b, the smaller the eps


```julia
prn{Tv,Ti}(G::SparseMatrixCSC{Tv,Ti}, s::Array{Int64,1}, phi::Float64, b::Int64)
```

localClustering.jl:366



### apr
computes an approximate page rank vector from a starting set s, an alpha and an epsilon   algorithm follows the Anderson,Chung,Lang paper and Dan Spielman's notes


```julia
apr{Tv,Ti}(G::SparseMatrixCSC{Tv,Ti}, s::Array{Int64,1}, alpha::Float64, eps::Float64)
```

localClustering.jl:432



