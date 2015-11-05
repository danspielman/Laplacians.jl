# cutPageRank
### pr
computes a page rank vector satisfying p = a/n * 1 + (1 - a) * W * p


```julia
pr{Tv,Ti}(G::SparseMatrixCSC{Tv,Ti}, alpha::Float64)
```

cutPageRank.jl:20



### ppr
computes the personal page rank vector from a starting vector s and an alpha; operates with lazy walk matrix 


```julia
ppr{Tv,Ti}(G::SparseMatrixCSC{Tv,Ti}, s::Array{Float64,1}, alpha::Float64)
ppr{Tv,Ti}(G::SparseMatrixCSC{Tv,Ti}, s::Array{Float64,1}, alpha::Float64, niter::Int64)
```

cutPageRank.jl:38



### apr
computes an approximate page rank vector from a starting vector s, an alpha and an epsilon 


```julia
apr{Tv,Ti}(G::SparseMatrixCSC{Tv,Ti}, s::Array{Float64,1}, alpha::Float64, eps::Float64)
```

cutPageRank.jl:71



### prn
prn{Tv, Ti}(G::SparseMatrixCSC{Tv, Ti}, v::Array{Int64,1}, phi::Float64, b::Int64)

the PageRank-Nibble cutting algorithm from the Anderson/Chung/Lang paper

v is a set of vertices, phi is a constant in (0, 1], and b is an integer in [1, [log m]]


```julia
prn{Tv,Ti}(G::SparseMatrixCSC{Tv,Ti}, v::Array{Int64,1}, phi::Float64, b::Int64)
```

cutPageRank.jl:115



