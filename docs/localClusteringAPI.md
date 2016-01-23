# localClustering
### localImprove
localImprove{Tv,Ti}(G::SparseMatrixCSC{Tv,Ti}, A::Array{Int64,1}; epsSigma=-1.0, err=1e-10, maxSize = max(G.n, G.m)

The LocalImprove function, from the Orrechia-Zhu paper. Given a graph and an initial set, finds a set of smaller conductance based on the starting set using a localized version of max-flow. 

G is the given graph, A is the initial set  epsSigma is a measure of the quality of the returning set (the smaller the better). It's defaulted to volume(A) / volume(VA) err is the numerical error considered throughout the algorithm. It's defaulted to 1e-10 maxSize is the maximum allowed size for the flow graph at any iteration of the algorithm. It's defaulted to |V|


```julia
localImprove{Tv,Ti}(G::SparseMatrixCSC{Tv,Ti}, A::Array{Int64,1})
```

localClustering.jl:17



### prn
prn{Tv, Ti}(G::SparseMatrixCSC{Tv,Ti}, s::Array{Int64,1}, phi::Float64, b::Int64)

The PageRank-Nibble cutting algorithm from the Anderson/Chung/Lang paper

s is a set of starting vertices, phi is a constant in (0, 1], and b is an integer in [1, [log m]]

phi is a bound on the quality of the conductance of the cut - the smaller the phi, the higher the quality.  b is used to handle precision throughout the algorithm - the higher the b, the greater the precision.


```julia
prn{Tv,Ti}(G::SparseMatrixCSC{Tv,Ti}, s::Array{Int64,1}, phi::Float64, b::Int64)
```

localClustering.jl:365



### apr
Computes an approximate page rank vector from a starting set s, an alpha and an epsilon The algorithm follows the Anderson,Chung,Lang paper and Dan Spielman's lecture notes


```julia
apr{Tv,Ti}(G::SparseMatrixCSC{Tv,Ti}, s::Array{Int64,1}, alpha::Float64, eps::Float64)
```

localClustering.jl:431



