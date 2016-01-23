### ErdosRenyi
Generate a random graph on n vertices with m edges. The actual number of edges will probably be smaller, as we sample with replacement


```julia
ErdosRenyi(n::Integer, m::Integer)
```

 at /Users/serbanstan/git/Laplacians.jl/src/graphGenerators.jl:337


### ErdosRenyiCluster
Generate an ER graph with average degree k, and then return the largest component. Will probably have fewer than n vertices. If you want to add a tree to bring it back to n, try ErdosRenyiClusterFix.


```julia
ErdosRenyiCluster(n::Integer, k::Integer)
```

 at /Users/serbanstan/git/Laplacians.jl/src/graphGenerators.jl:351


### ErdosRenyiClusterFix
Like an Erdos-Renyi cluster, but add back a tree so it has n vertices


```julia
ErdosRenyiClusterFix(n::Integer, k::Integer)
```

 at /Users/serbanstan/git/Laplacians.jl/src/graphGenerators.jl:364


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


### akpw
This is a wrapper for akpw!. It's slower, but won't modify the original graph. See akpw! documentation for more details.


```julia
akpw(origMat::SparseMatrixCSC{Tv,Ti<:Integer})
```

 at /Users/serbanstan/git/Laplacians.jl/src/akpwWeighted.jl:818


### akpw!
Constructs a low stretch tree using the Alon, Karp, Peleg, West algorithm. This version (akpw! instead of akpw) modifies the graph slightly changing the edges weights, then changing them back, which may lead to floating point imprecisions. akpw! is faster (about 10-20%), but akpw doesn't have float imprecisions.

The function has a few options:

kind: default is :max, which regards each edge weight as the inverse of its length (just like kruskal). If this is   set to anything else (e.g. :min), it will regard edge weight as length

randomClusters: default is false. This means the partition function searches for the beginning of the next cluster   in node order, rather than randomly choosing nodes. If this is set to false, it will   randomly choose the next node. This slows down akpw, but may produce better stretch.

metisClustering: default is false. If this is set to false, the graph will be partitioned   each time by metis, rather than by the akpw partitioning method.

shuffleClusters: default is true. This preserves the "reshuffleClusters" method after each each graph is   partitioned into clusters. If set to false, the function will skip this step. May be faster   but have worse stretch.

exponentialX: default is true, where the funciton exp(sqrt(log(nVertices) * log(log(nVertices)))) is used for X.   If set fo false, the function log(nVertices+1)/log(2) will be used for X instead. 

EXAMPLE:

[2, 1]  =  0.631273 [3, 1]  =  0.40103 [1, 2]  =  0.631273 [4, 2]  =  0.147018 [1, 3]  =  0.40103 [4, 3]  =  0.772661 [2, 4]  =  0.147018 [3, 4]  =  0.772661

```
  |
  |
  V
```

[2, 1]  =  0.631273 [3, 1]  =  0.40103 [1, 2]  =  0.631273 [1, 3]  =  0.40103 [4, 3]  =  0.772661 [3, 4]  =  0.772661


```julia
akpw!(mat::SparseMatrixCSC{Tv,Ti<:Integer})
```

 at /Users/serbanstan/git/Laplacians.jl/src/akpwWeighted.jl:733


### apr
Computes an approximate page rank vector from a starting set s, an alpha and an epsilon The algorithm follows the Anderson,Chung,Lang paper and Dan Spielman's lecture notes


```julia
apr{Tv,Ti}(G::SparseMatrixCSC{Tv,Ti}, s::Array{Int64,1}, alpha::Float64, eps::Float64)
```

 at /Users/serbanstan/git/Laplacians.jl/src/localClustering.jl:431


### augTreePrecon
This is an augmented spanning tree preconditioner for diagonally dominant linear systems.  It takes as optional input a tree growing algorithm. The default is a randomized variant of Kruskal. It adds back 2sqrt(n) edges via augmentTree. With the right tree, it should never be too bad.


```julia
augTreePrecon{Tv,Ti}(ddmat::SparseMatrixCSC{Tv,Ti})
```

 at /Users/serbanstan/git/Laplacians.jl/src/solvers.jl:188


### augTreeSolver
This is the solver that calls augTreePrecon


```julia
augTreeSolver{Tv,Ti}(ddmat::SparseMatrixCSC{Tv,Ti})
```

 at /Users/serbanstan/git/Laplacians.jl/src/solvers.jl:210


### augmentTree
Takes as input a tree and an adjacency matrix of a graph. It then computes the stretch of every edge of the graph wrt the tree.  It then adds back the k edges of highest stretch, and k edges sampled according to stretch


```julia
augmentTree{Tv,Ti}(tree::SparseMatrixCSC{Tv,Ti}, mat::SparseMatrixCSC{Tv,Ti}, k::Ti)
```

 at /Users/serbanstan/git/Laplacians.jl/src/solvers.jl:138


### backIndices
Same as the above, but now the graph is in adjacency list form 

Computes the back indices in a graph in O(M+N). works if for every edge (u,v), (v,u) is also in the graph 


```julia
backIndices{Tv,Ti}(G::SparseMatrixCSC{Tv,Ti})
backIndices{Tv1,Tv2}(G::Array{Array{Tuple{Tv1,Tv2},1},1})
```

 at /Users/serbanstan/git/Laplacians.jl/src/graphUtils.jl:35


### biggestComp
Return the biggest component in a graph, as a graph


```julia
biggestComp(mat::SparseMatrixCSC{Tv,Ti<:Integer})
```

 at /Users/serbanstan/git/Laplacians.jl/src/graphAlgs.jl:159


### cg

```julia
cg(mat, b::Array{Float64,1})
cg(mat, b::Array{Float32,1})
cg(mat, b)
```

 at /Users/serbanstan/git/Laplacians.jl/src/pcg.jl:29


### chimera
Builds the kth chimeric graph on n vertices. It does this by resetting the random number generator seed. It should captute the state of the generator before that and then return it, but it does not yet.

Builds a chimeric graph on n vertices. The components come from pureRandomGraph, connected by joinGraphs, productGraph and generalizedNecklace


```julia
chimera(n::Integer)
chimera(n::Integer, k::Integer)
```

 at /Users/serbanstan/git/Laplacians.jl/src/graphGenerators.jl:515


### compConductance
Returns the quality of the cut for a given graph and a given cut set s.   the result will be |outgoing edges| / min(|vertices in set|, |N - vertices in set|)


```julia
compConductance{Tv,Ti}(G::SparseMatrixCSC{Tv,Ti}, s::Array{Int64,1})
```

 at /Users/serbanstan/git/Laplacians.jl/src/graphUtils.jl:136


### compDepth

```julia
compDepth{Tv,Ti}(t::Laplacians.RootedTree{Tv,Ti})
```

 at /Users/serbanstan/git/Laplacians.jl/src/treeAlgs.jl:311


### compStretches
Compute the stretched of every edge in `mat` with respect to the tree `tree`. Returns the answer as a sparse matrix with the same nonzero structure as `mat`. Assumes that `mat` is symmetric. `tree` should be the adjacency matrix of a spanning tree.


```julia
compStretches{Tv,Ti}(t::Laplacians.RootedTree{Tv,Ti}, mat::SparseMatrixCSC{Tv,Ti})
compStretches{Tv,Ti}(tree::SparseMatrixCSC{Tv,Ti}, mat::SparseMatrixCSC{Tv,Ti})
```

 at /Users/serbanstan/git/Laplacians.jl/src/treeAlgs.jl:393


### completeBinaryTree
The complete binary tree on n vertices


```julia
completeBinaryTree(n::Int64)
```

 at /Users/serbanstan/git/Laplacians.jl/src/graphGenerators.jl:139


### completeGraph
The complete graph


```julia
completeGraph(n::Int64)
```

 at /Users/serbanstan/git/Laplacians.jl/src/graphGenerators.jl:17


### components
Computes the connected components of a graph. Returns them as a vector of length equal to the number of vertices. The vector numbers the components from 1 through the maximum number. For example,

```julia
gr = ErdosRenyi(10,11)
c = components(gr)

10-element Array{Int64,1}:
 1
 1
 1
 1
 2
 1
 1
 1
 3
 2
```


```julia
components{Tv,Ti}(mat::SparseMatrixCSC{Tv,Ti})
```

 at /Users/serbanstan/git/Laplacians.jl/src/graphAlgs.jl:65


### deg

```julia
deg{Tv,Ti}(mat::SparseMatrixCSC{Tv,Ti}, v::Ti)
```

 at /Users/serbanstan/git/Laplacians.jl/src/graphUtils.jl:11


### diagmat
Returns the diagonal matrix(as a sparse matrix) of a graph


```julia
diagmat{Tv,Ti}(G::SparseMatrixCSC{Tv,Ti})
```

 at /Users/serbanstan/git/Laplacians.jl/src/graphOps.jl:205


### dirEdgeVertexMat
The signed edge-vertex adjacency matrix


```julia
dirEdgeVertexMat(A::SparseMatrixCSC{Tv,Ti<:Integer})
```

 at /Users/serbanstan/git/Laplacians.jl/src/toposort.jl:49


### dumb
```
Modify a cluster by passing through all the vertices exactly once and 
adding/removing them based on the value of (Deg_external - Deg_Internal).
```


```julia
dumb{Tv,Ti}(G::SparseMatrixCSC{Tv,Ti}, s::Array{Int64,1})
```

 at /Users/serbanstan/git/Laplacians.jl/src/cutHeuristics.jl:106


### edgeVertexMat
The signed edge-vertex adjacency matrix


```julia
edgeVertexMat(mat::SparseMatrixCSC{Tv,Ti<:Integer})
```

 at /Users/serbanstan/git/Laplacians.jl/src/graphOps.jl:67


### findEntries
Similar to findnz, but also returns 0 entries that have an edge in the sparse matrix 


```julia
findEntries{Tv,Ti}(G::SparseMatrixCSC{Tv,Ti})
```

 at /Users/serbanstan/git/Laplacians.jl/src/graphUtils.jl:114


### floatGraph
Convert the nonzero entries in a graph to Float64


```julia
floatGraph(a::SparseMatrixCSC{Tv,Ti<:Integer})
```

 at /Users/serbanstan/git/Laplacians.jl/src/graphOps.jl:6


### generalizedNecklace
Constructs a generalized necklace graph starting with two graphs A and H. The resulting new graph will be constructed by expanding each vertex in H to an instance of A. k random edges will be generated between components. Thus, the resulting graph may have weighted edges.


```julia
generalizedNecklace{Tv,Ti}(A::SparseMatrixCSC{Tv,Ti}, H::SparseMatrixCSC{Tv,Ti<:Integer}, k::Int64)
```

 at /Users/serbanstan/git/Laplacians.jl/src/graphOps.jl:225


### generalizedRing
A generalization of a ring graph. The vertices are integers modulo n. Two are connected if their difference is in gens. For example, 

```
generalizedRing(17, [1 5])
```


```julia
generalizedRing(n::Int64, gens)
```

 at /Users/serbanstan/git/Laplacians.jl/src/graphGenerators.jl:38


### getObound
Computes the number of edges leaving s 


```julia
getObound{Tv,Ti}(G::SparseMatrixCSC{Tv,Ti}, s::Array{Int64,1})
```

 at /Users/serbanstan/git/Laplacians.jl/src/graphUtils.jl:167


### getVolume
Computes the volume of subset s in an unweighted graph G 


```julia
getVolume{Tv,Ti}(G::SparseMatrixCSC{Tv,Ti}, s::Array{Int64,1})
```

 at /Users/serbanstan/git/Laplacians.jl/src/graphUtils.jl:149


### grid2
An n-by-m grid graph.  iostropy is the weighting on edges in one direction.


```julia
grid2(n::Int64)
grid2(n::Int64, m::Int64)
```

 at /Users/serbanstan/git/Laplacians.jl/src/graphGenerators.jl:159


### grid2coords
Coordinates for plotting the vertices of the n-by-m grid graph


```julia
grid2coords(n::Int64, m::Int64)
grid2coords(n)
```

 at /Users/serbanstan/git/Laplacians.jl/src/graphGenerators.jl:163


### grownGraph
Create a graph on n vertices. For each vertex, give it k edges to randomly chosen prior vertices. This is a variety of a preferential attachment graph.    


```julia
grownGraph(n::Int64, k::Int64)
```

 at /Users/serbanstan/git/Laplacians.jl/src/graphGenerators.jl:202


### grownGraphD
Like a grownGraph, but it forces the edges to all be distinct. It starts out with a k+1 clique on the first k vertices


```julia
grownGraphD(n::Int64, k::Int64)
```

 at /Users/serbanstan/git/Laplacians.jl/src/graphGenerators.jl:234


### hyperCube
The d dimensional hypercube.  Has 2^d vertices


```julia
hyperCube(d::Int64)
```

 at /Users/serbanstan/git/Laplacians.jl/src/graphGenerators.jl:125


### isConnected
Returns true if graph is connected.  Calls components.


```julia
isConnected(mat::SparseMatrixCSC{Tv,Ti<:Integer})
```

 at /Users/serbanstan/git/Laplacians.jl/src/graphAlgs.jl:113


### joinGraphs
Create a disjoint union of graphs a and b,  and then put k random edges between them


```julia
joinGraphs{Tval,Tind}(a::SparseMatrixCSC{Tval,Tind}, b::SparseMatrixCSC{Tval,Tind}, k::Integer)
```

 at /Users/serbanstan/git/Laplacians.jl/src/graphOps.jl:111


### kruskal
Uses Kruskal's algorithm to compute a minimum (or maximum) spanning tree. Set kind=:max if you want the max spanning tree. It returns it a a graph


```julia
kruskal{Tv,Ti}(mat::SparseMatrixCSC{Tv,Ti})
```

 at /Users/serbanstan/git/Laplacians.jl/src/graphAlgs.jl:407


### lap
Create a Laplacian matrix from an adjacency matrix. We might want to do this differently, say by enforcing symmetry


```julia
lap(a)
```

 at /Users/serbanstan/git/Laplacians.jl/src/graphOps.jl:12


### lapChol


### lapWrapSolver
Takes a solver for solving nonsingular sdd systems, and returns a solver for solving Laplacian systems. The optional args tol and maxits are not necessarily taken by all solvers.  But, if they are, one can pass them here


```julia
lapWrapSolver(solver)
lapWrapSolver(solver, la::AbstractArray{T,N})
lapWrapSolver(solver, la::AbstractArray{T,N}, b)
```

 at /Users/serbanstan/git/Laplacians.jl/src/solvers.jl:108


### localImprove
localImprove{Tv,Ti}(G::SparseMatrixCSC{Tv,Ti}, A::Array{Int64,1}; epsSigma=-1.0, err=1e-10, maxSize = max(G.n, G.m)

The LocalImprove function, from the Orrechia-Zhu paper. Given a graph and an initial set, finds a set of smaller conductance based on the starting set using a localized version of max-flow. 

G is the given graph, A is the initial set  epsSigma is a measure of the quality of the returning set (the smaller the better). It's defaulted to volume(A) / volume(VA) err is the numerical error considered throughout the algorithm. It's defaulted to 1e-10 maxSize is the maximum allowed size for the flow graph at any iteration of the algorithm. It's defaulted to |V|


```julia
localImprove{Tv,Ti}(G::SparseMatrixCSC{Tv,Ti}, A::Array{Int64,1})
```

 at /Users/serbanstan/git/Laplacians.jl/src/localClustering.jl:17


### mapweight
Create a new graph that is the same as the original, but with f applied to each nonzero entry of a. For example, to make the weight of every edge uniform in [0,1], we could write

```julia
b = mapweight(a, x->rand(1)[1])
```


```julia
mapweight{Tval,Tind}(a::SparseMatrixCSC{Tval,Tind}, f)
```

 at /Users/serbanstan/git/Laplacians.jl/src/graphOps.jl:40


### matToTree

```julia
matToTree{Tv,Ti}(mat::SparseMatrixCSC{Tv,Ti})
matToTree{Tv,Ti}(mat::SparseMatrixCSC{Tv,Ti}, root::Ti)
```

 at /Users/serbanstan/git/Laplacians.jl/src/treeAlgs.jl:32


### matToTreeDepth

```julia
matToTreeDepth{Tv,Ti}(mat::SparseMatrixCSC{Tv,Ti})
matToTreeDepth{Tv,Ti}(mat::SparseMatrixCSC{Tv,Ti}, root::Ti)
```

 at /Users/serbanstan/git/Laplacians.jl/src/treeAlgs.jl:98


### maxflow
implementation of Dinic's algorithm. computes the maximum flow and min-cut in G between s and t    we consider the adjacency matrix to be the capacity matrix 


```julia
maxflow{Tv,Ti}(G::SparseMatrixCSC{Tv,Ti}, s::Int64, t::Int64)
```

 at /Users/serbanstan/git/Laplacians.jl/src/flow.jl:7


### nbri

```julia
nbri{Tv,Ti}(mat::SparseMatrixCSC{Tv,Ti}, v::Ti, i::Ti)
```

 at /Users/serbanstan/git/Laplacians.jl/src/graphUtils.jl:12


### nbrs

```julia
nbrs{Tv,Ti}(mat::SparseMatrixCSC{Tv,Ti}, v::Ti)
```

 at /Users/serbanstan/git/Laplacians.jl/src/graphUtils.jl:14


### pathFromParents

```julia
pathFromParents(parents, y)
```

 at /Users/serbanstan/git/Laplacians.jl/src/graphAlgs.jl:215


### pathGraph
The path graph on n vertices


```julia
pathGraph(n::Int64)
```

 at /Users/serbanstan/git/Laplacians.jl/src/graphGenerators.jl:8


### pcg

```julia
pcg(mat, b::Array{Float64,1}, pre)
pcg(mat, b::Array{Float32,1}, pre)
pcg(mat, b, pre)
```

 at /Users/serbanstan/git/Laplacians.jl/src/pcg.jl:42


### plotGraph
Plots graph gr with coordinates (x,y)


```julia
plotGraph(gr, x, y)
plotGraph(gr, x, y, color)
```

 at /Users/serbanstan/git/Laplacians.jl/src/graphOps.jl:129


### prefAttach
A preferential attachment graph in which each vertex has k edges to those that come before.  These are chosen with probability p to be from a random vertex, and with probability 1-p to come from the endpoint of a random edge. It begins with a k-clique on the first k+1 vertices.


```julia
prefAttach(n::Int64, k::Int64, p::Float64)
```

 at /Users/serbanstan/git/Laplacians.jl/src/graphGenerators.jl:258


### prim

```julia
prim(mat::SparseMatrixCSC{Tv,Ti<:Integer})
```

 at /Users/serbanstan/git/Laplacians.jl/src/graphAlgs.jl:438


### prn
prn{Tv, Ti}(G::SparseMatrixCSC{Tv,Ti}, s::Array{Int64,1}, phi::Float64, b::Int64)

The PageRank-Nibble cutting algorithm from the Anderson/Chung/Lang paper

s is a set of starting vertices, phi is a constant in (0, 1], and b is an integer in [1, [log m]]

phi is a bound on the quality of the conductance of the cut - the smaller the phi, the higher the quality.  b is used to handle precision throughout the algorithm - the higher the b, the greater the precision.


```julia
prn{Tv,Ti}(G::SparseMatrixCSC{Tv,Ti}, s::Array{Int64,1}, phi::Float64, b::Int64)
```

 at /Users/serbanstan/git/Laplacians.jl/src/localClustering.jl:365


### productGraph
The Cartesian product of two graphs.  When applied to two paths, it gives a grid.


```julia
productGraph(a0::SparseMatrixCSC{Tv,Ti<:Integer}, a1::SparseMatrixCSC{Tv,Ti<:Integer})
```

 at /Users/serbanstan/git/Laplacians.jl/src/graphOps.jl:58


### pureRandomGraph
Generate a random graph with n vertices from one of our natural distributions


```julia
pureRandomGraph(n::Integer)
```

 at /Users/serbanstan/git/Laplacians.jl/src/graphGenerators.jl:379


### randGenRing
A random generalized ring graph of degree k. Gens always contains 1, and the other k-1 edge types are chosen from an exponential distribution


```julia
randGenRing(n::Int64, k::Integer)
```

 at /Users/serbanstan/git/Laplacians.jl/src/graphGenerators.jl:62


### randMatching
A random matching on n vertices


```julia
randMatching(n::Int64)
```

 at /Users/serbanstan/git/Laplacians.jl/src/graphGenerators.jl:174


### randRegular
A sum of k random matchings on n vertices


```julia
randRegular(n::Int64, k::Int64)
```

 at /Users/serbanstan/git/Laplacians.jl/src/graphGenerators.jl:187


### randWeight
Applies one of a number of random weighting schemes to the edges of the graph


```julia
randWeight(a)
```

 at /Users/serbanstan/git/Laplacians.jl/src/graphGenerators.jl:535


### randishKruskal

```julia
randishKruskal{Tv,Ti}(mat::SparseMatrixCSC{Tv,Ti})
```

 at /Users/serbanstan/git/Laplacians.jl/src/randTrees.jl:10


### randishPrim

```julia
randishPrim{Tval,Tind}(mat::SparseMatrixCSC{Tval,Tind})
```

 at /Users/serbanstan/git/Laplacians.jl/src/randTrees.jl:47


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
To read a simple edge list, each line being an (i, j) pair


```julia
readIJ(filename::AbstractString)
readIJ(filename::AbstractString, sep)
```

 at /Users/serbanstan/git/Laplacians.jl/src/IO.jl:4


### readIJV
To read a simple edge list, each line being an (i, j, v) pair. The parens should not be there in the format, just commas separating. To generate this format in Matlab, you just need to be careful to write the vertex indices with sufficient precision.  For example, you can do this

```
>> [ai,aj,av] = find(triu(a));
>> dlmwrite('graph.txt',[ai,aj,av],'precision',9);
```


```julia
readIJV(filename::AbstractString)
```

 at /Users/serbanstan/git/Laplacians.jl/src/IO.jl:25


### refineCut
```
Modifies a cluster by adding or removing vertices by picking at each step 
the vertex that has the maximum value of (Deg_external - Deg_Internal).
Each vertex can be added in/removed only once.
```


```julia
refineCut{Tv,Ti}(G::SparseMatrixCSC{Tv,Ti}, s::Array{Int64,1})
```

 at /Users/serbanstan/git/Laplacians.jl/src/cutHeuristics.jl:9


### ringGraph
The simple ring on n vertices


```julia
ringGraph(n::Int64)
```

 at /Users/serbanstan/git/Laplacians.jl/src/graphGenerators.jl:23


### semiWtedChimera
A Chimera graph with some weights.  The weights just appear when graphs are combined. For more interesting weights, use `wtedChimera`


```julia
semiWtedChimera(n::Integer)
semiWtedChimera(n::Integer, k::Integer)
```

 at /Users/serbanstan/git/Laplacians.jl/src/graphGenerators.jl:442


### setValue
Sets the value of a certain edge in a sparse graph; value can be 0 without the edges dissapearing 


```julia
setValue{Tv,Ti}(mat::SparseMatrixCSC{Tv,Ti}, v::Ti, i::Ti, a::Tv)
```

 at /Users/serbanstan/git/Laplacians.jl/src/graphUtils.jl:29


### shortIntGraph
Convert the indices in a graph to 32-bit ints.  This takes less storage, but does not speed up much


```julia
shortIntGraph(a::SparseMatrixCSC{Tv,Ti<:Integer})
```

 at /Users/serbanstan/git/Laplacians.jl/src/graphOps.jl:3


### shortestPathTree
Computes the shortest path tree, and returns it as a sparse matrix. Treats edge weights as reciprocals of lengths. For example:

```julia
a = [0 2 1; 2 0 3; 1 3 0]
tr = full(shortestPathTree(sparse(a),1))

3x3 Array{Float64,2}:
 0.0  2.0  0.0
 2.0  0.0  3.0
 0.0  3.0  0.0
```


```julia
shortestPathTree(a, start)
```

 at /Users/serbanstan/git/Laplacians.jl/src/graphAlgs.jl:239


### shortestPaths
Computes the lenghts of shortest paths from `start`. Returns both a vector of the lenghts, and the parent array in the shortest path tree.

This algorithm treats edge weights as reciprocals of distances. DOC BETTER


```julia
shortestPaths{Tv,Ti}(mat::SparseMatrixCSC{Tv,Ti}, start::Ti)
shortestPaths{Tv,Ti}(mat::SparseMatrixCSC{Tv,Ti})
```

 at /Users/serbanstan/git/Laplacians.jl/src/graphAlgs.jl:175


### spectralCoords
Computes the spectral coordinates of a graph


```julia
spectralCoords(a)
```

 at /Users/serbanstan/git/Laplacians.jl/src/graphOps.jl:174


### spectralDrawing
Computes spectral coordinates, and then uses plotGraph to draw


```julia
spectralDrawing(a)
```

 at /Users/serbanstan/git/Laplacians.jl/src/graphOps.jl:166


### subsampleEdges
Create a new graph from the old, but keeping edge edge with probability `p`


```julia
subsampleEdges(a::SparseMatrixCSC{Float64,Int64}, p::Float64)
```

 at /Users/serbanstan/git/Laplacians.jl/src/graphOps.jl:76


### tarjanStretch

```julia
tarjanStretch{Tv,Ti}(t::Laplacians.RootedTree{Tv,Ti}, mat::SparseMatrixCSC{Tv,Ti}, depth::Array{Tv,1})
```

 at /Users/serbanstan/git/Laplacians.jl/src/treeAlgs.jl:334


### toUnitVector
Creates a unit vector of length n from a given set of integers, with weights based on the number of occurences


```julia
toUnitVector(a::Array{Int64,1}, n)
```

 at /Users/serbanstan/git/Laplacians.jl/src/graphOps.jl:183


### toposort

```julia
toposort{Tv,Ti}(mat::SparseMatrixCSC{Tv,Ti})
```

 at /Users/serbanstan/git/Laplacians.jl/src/toposort.jl:13


### twoLift
Creats a 2-lift of a.  `flip` is a boolean indicating which edges cross


```julia
twoLift(a)
twoLift(a, flip::AbstractArray{Bool,1})
twoLift(a, k::Integer)
```

 at /Users/serbanstan/git/Laplacians.jl/src/graphOps.jl:99


### uniformWeight
Put a uniform [0,1] weight on every edge.  This is an example of how to use mapweight.


```julia
uniformWeight{Tval,Tind}(a::SparseMatrixCSC{Tval,Tind})
```

 at /Users/serbanstan/git/Laplacians.jl/src/graphOps.jl:49


### uniformWeight!
Set the weight of every edge to 1


```julia
uniformWeight!(mat::SparseMatrixCSC{Tv,Ti<:Integer})
```

 at /Users/serbanstan/git/Laplacians.jl/src/graphOps.jl:53


### unweight
Create a new graph in that is the same as the original, but with all edge weights 1


```julia
unweight{Tval,Tind}(ain::SparseMatrixCSC{Tval,Tind})
```

 at /Users/serbanstan/git/Laplacians.jl/src/graphOps.jl:16


### unweight!
Change the weight of every edge in a to 1


```julia
unweight!{Tval,Tind}(a::SparseMatrixCSC{Tval,Tind})
```

 at /Users/serbanstan/git/Laplacians.jl/src/graphOps.jl:26


### vecToComps
This turns a component vector, like that generated by components, into an array of arrays of indices of vertices in each component.  For example,

```julia
comps = vecToComps(c)

3-element Array{Array{Int64,1},1}:
 [1,2,3,4,6,7,8]
 [5,10]
 [9]
```


```julia
vecToComps{Ti}(compvec::Array{Ti,1})
```

 at /Users/serbanstan/git/Laplacians.jl/src/graphAlgs.jl:136


### wdeg
Finds the weighted degree of a vertex in the graph 


```julia
wdeg{Tv,Ti}(mat::SparseMatrixCSC{Tv,Ti}, v::Ti)
```

 at /Users/serbanstan/git/Laplacians.jl/src/graphUtils.jl:19


### weighti

```julia
weighti{Tv,Ti}(mat::SparseMatrixCSC{Tv,Ti}, v::Ti, i::Ti)
```

 at /Users/serbanstan/git/Laplacians.jl/src/graphUtils.jl:13


### writeIJV
Writes the upper portion of a matrix in ijv format, one row for each edge, separated by commas.  Only writes the upper triangular portion. The result can be read from Matlab like this:

```
>> dl = dlmread('graph.txt');
>> a = sparse(dl(:,1),dl(:,2),dl(:,3));
>> n = max(size(a))
>> a(n,n) = 0;
>> a = a + a';
```


```julia
writeIJV(filename::AbstractString, mat)
```

 at /Users/serbanstan/git/Laplacians.jl/src/IO.jl:52


### wtedChimera
Builds the kth wted chimeric graph on n vertices. It does this by resetting the random number generator seed. It should captute the state of the generator before that and then return it, but it does not yet.

Generate a chimera, and then apply a random weighting scheme


```julia
wtedChimera(n::Integer)
wtedChimera(n::Integer, k::Integer)
```

 at /Users/serbanstan/git/Laplacians.jl/src/graphGenerators.jl:607


