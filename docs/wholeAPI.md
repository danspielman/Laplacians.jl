### ErdosRenyiCluster
Generate an ER graph with average degree k, and then return the largest component. Will probably have fewer than n vertices. If you want to add a tree to bring it back to n, try ErdosRenyiClusterFix.


```julia
ErdosRenyiCluster(n::Integer, k::Integer)
```

 at /Users/spielman/git/Laplacians.jl/src/graphGenerators.jl:319


### ErdosRenyiClusterFix
Like an Erdos-Renyi cluster, but add back a tree so it has n vertices


```julia
ErdosRenyiClusterFix(n::Integer, k::Integer)
```

 at /Users/spielman/git/Laplacians.jl/src/graphGenerators.jl:332


### Laplacians
A package for graph computations related to graph Laplacians

Graphs are represented by sparse adjacency matrices, etc.



### RootedTree
**Summary:**

```julia
type Laplacians.RootedTree{Tval,Tind} <: Any
```

**Fields:**

```julia
root     :: Tind
parent   :: Array{Tind,1}
children :: Array{Tind,1}
weights  :: Array{Tval,1}
numKids  :: Array{Tind,1}
kidsPtr  :: Array{Tind,1}
```


### apr
computes an approximate page rank vector from a starting vector s, an alpha and an epsilon 


```julia
apr{Tv,Ti}(G::SparseMatrixCSC{Tv,Ti}, s::Array{Float64,1}, alpha::Float64, eps::Float64)
```

 at /Users/spielman/git/Laplacians.jl/src/cutPageRank.jl:71


### backIndices
computes the back indices in a graph in O(M+N) 


```julia
backIndices{Tv,Ti}(G::SparseMatrixCSC{Tv,Ti})
```

 at /Users/spielman/git/Laplacians.jl/src/graphUtils.jl:22


### biggestComp
Return the biggest component in a graph


```julia
biggestComp(mat::SparseMatrixCSC{Tv,Ti<:Integer})
```

 at /Users/spielman/git/Laplacians.jl/src/graphAlgs.jl:111


### cg

```julia
cg(mat, b::Array{Float64,1})
cg(mat, b::Array{Float32,1})
cg(mat, b)
```

 at /Users/spielman/git/Laplacians.jl/src/pcg.jl:29


### chimera
Builds the kth chimeric graph on n vertices. It does this by resetting the random number generator seed. It should captute the state of the generator before that and then return it, but it does not yet.

Builds a chimeric graph on n vertices. The components come from pureRandomGraph, connected by joinGraphs, productGraph and generalizedNecklace


```julia
chimera(n::Integer)
chimera(n::Integer, k::Integer)
```

 at /Users/spielman/git/Laplacians.jl/src/graphGenerators.jl:392


### compDepth

```julia
compDepth{Tv,Ti}(t::Laplacians.RootedTree{Tv,Ti})
```

 at /Users/spielman/git/Laplacians.jl/src/treeAlgs.jl:249


### compStretches

```julia
compStretches{Tv,Ti}(t::Laplacians.RootedTree{Tv,Ti}, mat::SparseMatrixCSC{Tv,Ti})
compStretches{Tv,Ti}(tree::SparseMatrixCSC{Tv,Ti}, mat::SparseMatrixCSC{Tv,Ti})
```

 at /Users/spielman/git/Laplacians.jl/src/treeAlgs.jl:319


### completeBinaryTree
The complete binary tree on n vertices


```julia
completeBinaryTree(n::Int64)
```

 at /Users/spielman/git/Laplacians.jl/src/graphGenerators.jl:119


### completeGraph
The complete graph


```julia
completeGraph(n::Int64)
```

 at /Users/spielman/git/Laplacians.jl/src/graphGenerators.jl:54


### components

```julia
components{Tv,Ti}(mat::SparseMatrixCSC{Tv,Ti})
```

 at /Users/spielman/git/Laplacians.jl/src/graphAlgs.jl:43


### deg

```julia
deg{Tv,Ti}(mat::SparseMatrixCSC{Tv,Ti}, v::Ti)
```

 at /Users/spielman/git/Laplacians.jl/src/graphUtils.jl:11


### diagmat
returns the diagonal matrix(as a sparse matrix) of a graph


```julia
diagmat{Tv,Ti}(G::SparseMatrixCSC{Tv,Ti})
```

 at /Users/spielman/git/Laplacians.jl/src/graphOps.jl:214


### edgeVertexMat
The signed edge-vertex adjacency matrix


```julia
edgeVertexMat(mat::SparseMatrixCSC{Tv,Ti<:Integer})
```

 at /Users/spielman/git/Laplacians.jl/src/graphOps.jl:78


### generalizedNecklace
Constructs a generalized necklace graph starting with two graphs A and H. The resulting new graph will be constructed by expanding each vertex in H to an instance of A. k random edges will be generated between components. Thus, the resulting graph may have weighted edges.


```julia
generalizedNecklace{Tv,Ti}(A::SparseMatrixCSC{Tv,Ti}, H::SparseMatrixCSC{Tv,Ti<:Integer}, k::Int64)
```

 at /Users/spielman/git/Laplacians.jl/src/graphOps.jl:234


### generalizedRing
A generalization of a ring graph. The vertices are integers modulo n. Two are connected if their difference is in gens. For example, 

```
generalizedRing(17, [1 5])
```


```julia
generalizedRing(n::Int64, gens)
```

 at /Users/spielman/git/Laplacians.jl/src/graphGenerators.jl:75


### grid2
An n-by-m grid graph.  iostropy is the weighting on edges in one direction.


```julia
grid2(n::Int64)
grid2(n::Int64, m::Int64)
```

 at /Users/spielman/git/Laplacians.jl/src/graphGenerators.jl:139


### grid2coords
Coordinates for plotting the vertices of the n-by-m grid graph


```julia
grid2coords(n::Int64, m::Int64)
grid2coords(n)
```

 at /Users/spielman/git/Laplacians.jl/src/graphGenerators.jl:143


### grownGraph
Create a graph on n vertices. For each vertex, give it k edges to randomly chosen prior vertices. This is a variety of a preferential attachment graph.    


```julia
grownGraph(n::Int64, k::Int64)
```

 at /Users/spielman/git/Laplacians.jl/src/graphGenerators.jl:182


### grownGraphD
Like a grownGraph, but it forces the edges to all be distinct. It starts out with a k+1 clique on the first k vertices


```julia
grownGraphD(n::Int64, k::Int64)
```

 at /Users/spielman/git/Laplacians.jl/src/graphGenerators.jl:214


### hyperCube
The d dimensional hypercube.  Has 2^d vertices


```julia
hyperCube(d::Int64)
```

 at /Users/spielman/git/Laplacians.jl/src/graphGenerators.jl:105


### joinGraphs
create a disjoint union of graphs a and b,  and then put k random edges between them


```julia
joinGraphs{Tval,Tind}(a::SparseMatrixCSC{Tval,Tind}, b::SparseMatrixCSC{Tval,Tind}, k::Integer)
```

 at /Users/spielman/git/Laplacians.jl/src/graphOps.jl:120


### kruskal

```julia
kruskal{Tv,Ti}(mat::SparseMatrixCSC{Tv,Ti})
```

 at /Users/spielman/git/Laplacians.jl/src/graphAlgs.jl:316


### lap
Create a Laplacian matrix from an adjacency matrix. We might want to do this differently, say by enforcing symmetry


```julia
lap(a)
```

 at /Users/spielman/git/Laplacians.jl/src/graphOps.jl:39


### mapweight
Create a new graph that is the same as the original, but with f applied to each nonzero entry of a. For example, to make the weight of every edge uniform in [0,1], we could write

```julia
b = mapweight(a, x->rand(1)[1])
```


```julia
mapweight{Tval,Tind}(a::SparseMatrixCSC{Tval,Tind}, f)
```

 at /Users/spielman/git/Laplacians.jl/src/graphOps.jl:56


### matToTree

```julia
matToTree{Tv,Ti}(mat::SparseMatrixCSC{Tv,Ti})
matToTree{Tv,Ti}(mat::SparseMatrixCSC{Tv,Ti}, root::Ti)
```

 at /Users/spielman/git/Laplacians.jl/src/treeAlgs.jl:32


### matToTreeDepth

```julia
matToTreeDepth{Tv,Ti}(mat::SparseMatrixCSC{Tv,Ti})
matToTreeDepth{Tv,Ti}(mat::SparseMatrixCSC{Tv,Ti}, root::Ti)
```

 at /Users/spielman/git/Laplacians.jl/src/treeAlgs.jl:98


### maxflow
implementation of Dinic's algorithm. computes the maximum flow and min-cut in G between s and t. we consider the adjacency matrix to be the capacity matrix 


```julia
maxflow{Tv,Ti}(G::SparseMatrixCSC{Tv,Ti}, s::Int64, t::Int64)
```

 at /Users/spielman/git/Laplacians.jl/src/flow.jl:4


### nbri

```julia
nbri{Tv,Ti}(mat::SparseMatrixCSC{Tv,Ti}, v::Ti, i::Ti)
```

 at /Users/spielman/git/Laplacians.jl/src/graphUtils.jl:12


### nbrs

```julia
nbrs{Tv,Ti}(mat::SparseMatrixCSC{Tv,Ti}, v::Ti)
```

 at /Users/spielman/git/Laplacians.jl/src/graphUtils.jl:14


### pathGraph
The path graph on n vertices


```julia
pathGraph(n::Int64)
```

 at /Users/spielman/git/Laplacians.jl/src/graphGenerators.jl:45


### pcg

```julia
pcg(mat, b::Array{Float64,1}, pre)
pcg(mat, b::Array{Float32,1}, pre)
pcg(mat, b, pre)
```

 at /Users/spielman/git/Laplacians.jl/src/pcg.jl:42


### plotGraph
Plots graph gr with coordinates (x,y)


```julia
plotGraph(gr, x, y)
plotGraph(gr, x, y, color)
```

 at /Users/spielman/git/Laplacians.jl/src/graphOps.jl:138


### ppr
computes the personal page rank vector from a starting vector s and an alpha; operates with lazy walk matrix 


```julia
ppr{Tv,Ti}(G::SparseMatrixCSC{Tv,Ti}, s::Array{Float64,1}, alpha::Float64)
ppr{Tv,Ti}(G::SparseMatrixCSC{Tv,Ti}, s::Array{Float64,1}, alpha::Float64, niter::Int64)
```

 at /Users/spielman/git/Laplacians.jl/src/cutPageRank.jl:38


### pr
computes a page rank vector satisfying p = a/n * 1 + (1 - a) * W * p


```julia
pr{Tv,Ti}(G::SparseMatrixCSC{Tv,Ti}, alpha::Float64)
```

 at /Users/spielman/git/Laplacians.jl/src/cutPageRank.jl:20


### prefAttach
A preferential attachment graph in which each vertex has k edges to those that come before.  These are chosen with probability p to be from a random vertex, and with probability 1-p to come from the endpoint of a random edge. It begins with a k-clique on the first k+1 vertices.


```julia
prefAttach(n::Int64, k::Int64, p::Float64)
```

 at /Users/spielman/git/Laplacians.jl/src/graphGenerators.jl:238


### prn
prn{Tv, Ti}(G::SparseMatrixCSC{Tv, Ti}, v::Array{Int64,1}, phi::Float64, b::Int64)

the PageRank-Nibble cutting algorithm from the Anderson/Chung/Lang paper

v is a set of vertices, phi is a constant in (0, 1], and b is an integer in [1, [log m]]


```julia
prn{Tv,Ti}(G::SparseMatrixCSC{Tv,Ti}, v::Array{Int64,1}, phi::Float64, b::Int64)
```

 at /Users/spielman/git/Laplacians.jl/src/cutPageRank.jl:115


### productGraph
The Cartesian product of two graphs.  When applied to two paths, it gives a grid.


```julia
productGraph(a0::SparseMatrixCSC{Tv,Ti<:Integer}, a1::SparseMatrixCSC{Tv,Ti<:Integer})
```

 at /Users/spielman/git/Laplacians.jl/src/graphOps.jl:69


### pureRandomGraph
Generate a random graph with n vertices from one of our natural distributions


```julia
pureRandomGraph(n::Integer)
```

 at /Users/spielman/git/Laplacians.jl/src/graphGenerators.jl:347


### randGenRing
A random generalized ring graph of degree k. Gens always contains 1, and the other k-1 edge types are chosen from an exponential distribution


```julia
randGenRing(n::Int64, k::Integer)
```

 at /Users/spielman/git/Laplacians.jl/src/graphGenerators.jl:98


### randMatching
A random matching on n vertices


```julia
randMatching(n::Int64)
```

 at /Users/spielman/git/Laplacians.jl/src/graphGenerators.jl:154


### randRegular
A sum of k random matchings on n vertices


```julia
randRegular(n::Int64, k::Int64)
```

 at /Users/spielman/git/Laplacians.jl/src/graphGenerators.jl:167


### randWeight
Applies one of a number of random weighting schemes to the edges of the graph


```julia
randWeight(a)
```

 at /Users/spielman/git/Laplacians.jl/src/graphGenerators.jl:463


### randperm
```rst
..  randperm([rng,] n)

Construct a random permutation of length ``n``. The optional ``rng`` argument
specifies a random number generator, see :ref:`Random Numbers <random-numbers>`.
```

Randomly permutes the vertex indices


```julia
randperm(r::AbstractRNG, n::Integer)
randperm(n::Integer)
randperm(mat::AbstractArray{T,2})
randperm(f::Expr)
```

 at random.jl:1341


### readIJ

```julia
readIJ(filename::AbstractString)
```

 at /Users/spielman/git/Laplacians.jl/src/graphGenerators.jl:7


### ringGraph
The simple ring on n vertices


```julia
ringGraph(n::Int64)
```

 at /Users/spielman/git/Laplacians.jl/src/graphGenerators.jl:60


### setValue

```julia
setValue{Tv,Ti}(mat::SparseMatrixCSC{Tv,Ti}, v::Ti, i::Ti, a::Tv)
```

 at /Users/spielman/git/Laplacians.jl/src/graphUtils.jl:17


### shortIntGraph
Convert the indices in a graph to 32-bit ints.  This takes less storage, but does not speed up much


```julia
shortIntGraph(a::SparseMatrixCSC{Tv,Ti<:Integer})
```

 at /Users/spielman/git/Laplacians.jl/src/graphOps.jl:30


### shortestPaths

```julia
shortestPaths{Tv,Ti}(mat::SparseMatrixCSC{Tv,Ti}, start::Ti)
```

 at /Users/spielman/git/Laplacians.jl/src/graphAlgs.jl:122


### spectralCoords
Computes the spectral coordinates of a graph


```julia
spectralCoords(a)
```

 at /Users/spielman/git/Laplacians.jl/src/graphOps.jl:183


### spectralDrawing
Computes spectral coordinates, and then uses plotGraph to draw


```julia
spectralDrawing(a)
```

 at /Users/spielman/git/Laplacians.jl/src/graphOps.jl:175


### subsampleEdges
Create a new graph from the old, but keeping edge edge with probability `p`


```julia
subsampleEdges(a::SparseMatrixCSC{Float64,Int64}, p::Float64)
```

 at /Users/spielman/git/Laplacians.jl/src/graphOps.jl:87


### tarjanStretch

```julia
tarjanStretch{Tv,Ti}(t::Laplacians.RootedTree{Tv,Ti}, mat::SparseMatrixCSC{Tv,Ti}, depth::Array{Tv,1})
```

 at /Users/spielman/git/Laplacians.jl/src/treeAlgs.jl:272


### toUnitVector
creates a unit vector of length n from a given set of integers, with weights based on the number of occurences


```julia
toUnitVector(a::Array{Int64,1}, n)
```

 at /Users/spielman/git/Laplacians.jl/src/graphOps.jl:192


### twoLift
Creats a 2-lift of a.  `flip` is a boolean indicating which edges cross


```julia
twoLift(a)
twoLift(a, flip::AbstractArray{Bool,1})
twoLift(a, k::Integer)
```

 at /Users/spielman/git/Laplacians.jl/src/graphOps.jl:108


### uniformWeight
Put a uniform [0,1] weight on every edge.  This is an example of how to use mapweight.


```julia
uniformWeight{Tval,Tind}(a::SparseMatrixCSC{Tval,Tind})
```

 at /Users/spielman/git/Laplacians.jl/src/graphOps.jl:65


### unweight
Create a new graph in that is the same as the original, but with all edge weights 1


```julia
unweight{Tval,Tind}(a::SparseMatrixCSC{Tval,Tind})
```

 at /Users/spielman/git/Laplacians.jl/src/graphOps.jl:43


### vecToComps

```julia
vecToComps{Ti}(compvec::Array{Ti,1})
```

 at /Users/spielman/git/Laplacians.jl/src/graphAlgs.jl:88


### weighti

```julia
weighti{Tv,Ti}(mat::SparseMatrixCSC{Tv,Ti}, v::Ti, i::Ti)
```

 at /Users/spielman/git/Laplacians.jl/src/graphUtils.jl:13


### wtedChimera
Builds the kth wted chimeric graph on n vertices. It does this by resetting the random number generator seed. It should captute the state of the generator before that and then return it, but it does not yet.

Generate a chimera, and then apply a random weighting scheme


```julia
wtedChimera(n::Integer)
wtedChimera(n::Integer, k::Integer)
```

 at /Users/spielman/git/Laplacians.jl/src/graphGenerators.jl:528


