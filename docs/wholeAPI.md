### ErdosRenyi
Generate a random graph on n vertices with m edges. The actual number of edges will probably be smaller, as we sample with replacement


```julia
ErdosRenyi(n::Integer, m::Integer)
```

 at /Users/spielman/.julia/v0.4/Laplacians/src/graphGenerators.jl:337


### ErdosRenyiCluster
Generate an ER graph with average degree k, and then return the largest component. Will probably have fewer than n vertices. If you want to add a tree to bring it back to n, try ErdosRenyiClusterFix.


```julia
ErdosRenyiCluster(n::Integer, k::Integer)
```

 at /Users/spielman/.julia/v0.4/Laplacians/src/graphGenerators.jl:351


### ErdosRenyiClusterFix
Like an Erdos-Renyi cluster, but add back a tree so it has n vertices


```julia
ErdosRenyiClusterFix(n::Integer, k::Integer)
```

 at /Users/spielman/.julia/v0.4/Laplacians/src/graphGenerators.jl:364


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
Computes a low stretch spanning tree of `graph`, and returns it as a graph. The default version is 0.  In event of emergency, one can try `ver=2`.  It is usually slower, but might have slightly better stretch.


```julia
akpw(graph)
```

 at /Users/spielman/.julia/v0.4/Laplacians/src/akpw.jl:316


### akpwU
Computes a low stretch spanning tree of an unweighted `graph`, and returns it as a graph.


```julia
akpwU(graph)
```

 at /Users/spielman/.julia/v0.4/Laplacians/src/akpw.jl:166


### apr
Computes an approximate page rank vector from a starting set s, an alpha and an epsilon The algorithm follows the Anderson,Chung,Lang paper and Dan Spielman's lecture notes


```julia
apr{Tv,Ti}(G::SparseMatrixCSC{Tv,Ti}, s::Array{Int64,1}, alpha::Float64, eps::Float64)
```

 at /Users/spielman/.julia/v0.4/Laplacians/src/localClustering.jl:440


### augTreeLapPrecon
A version of augTreePrecon specialized for Laplacians


```julia
augTreeLapPrecon{Tv,Ti}(ddmat::SparseMatrixCSC{Tv,Ti})
```

 at /Users/spielman/.julia/v0.4/Laplacians/src/solvers.jl:291


### augTreeLapSolver
This is the solver that calls augTreeLapPrecon.  A solver specialized for Laplacians.


```julia
augTreeLapSolver{Tv,Ti}(ddmat::SparseMatrixCSC{Tv,Ti})
```

 at /Users/spielman/.julia/v0.4/Laplacians/src/solvers.jl:313


### augTreePrecon
This is an augmented spanning tree preconditioner for diagonally dominant linear systems.  It takes as optional input a tree growing algorithm. The default is a randomized variant of Kruskal. It adds back 2sqrt(n) edges via augmentTree. With the right tree, it should never be too bad.


```julia
augTreePrecon{Tv,Ti}(ddmat::SparseMatrixCSC{Tv,Ti})
```

 at /Users/spielman/.julia/v0.4/Laplacians/src/solvers.jl:258


### augTreeSolver
This is the solver that calls augTreePrecon


```julia
augTreeSolver{Tv,Ti}(ddmat::SparseMatrixCSC{Tv,Ti})
```

 at /Users/spielman/.julia/v0.4/Laplacians/src/solvers.jl:280


### augmentTree
Takes as input a tree and an adjacency matrix of a graph. It then computes the stretch of every edge of the graph wrt the tree.  It then adds back the k edges of highest stretch, and k edges sampled according to stretch


```julia
augmentTree{Tv,Ti}(tree::SparseMatrixCSC{Tv,Ti}, mat::SparseMatrixCSC{Tv,Ti}, k::Ti)
```

 at /Users/spielman/.julia/v0.4/Laplacians/src/solvers.jl:208


### backIndices
Same as the above, but now the graph is in adjacency list form 

Computes the back indices in a graph in O(M+N). works if for every edge (u,v), (v,u) is also in the graph 


```julia
backIndices{Tv,Ti}(G::SparseMatrixCSC{Tv,Ti})
backIndices{Tv1,Tv2}(G::Array{Array{Tuple{Tv1,Tv2},1},1})
```

 at /Users/spielman/.julia/v0.4/Laplacians/src/graphUtils.jl:62


### biggestComp
Return the biggest component in a graph, as a graph


```julia
biggestComp(mat::SparseMatrixCSC{Tv,Ti<:Integer})
```

 at /Users/spielman/.julia/v0.4/Laplacians/src/graphAlgs.jl:159


### cg
`cg(mat, b; tol, maxits, maxtime, verbose)` solves a symmetric linear system. `tol` is set to 1e-6 by default, `maxits` defaults to Inf `maxtime` defaults to Inf.  It measures seconds. `verbose` defaults to false


```julia
cg(mat, b::Array{Float64,1})
cg(mat, b::Array{Float32,1})
cg(mat, b)
```

 at /Users/spielman/.julia/v0.4/Laplacians/src/pcg.jl:45


### cgSolver
`cgSolver(mat; tol, maxits, maxtime, verbose)` creates a solver for a PSD system `mat`. The parameters are as described in cg.


```julia
cgSolver(mat)
```

 at /Users/spielman/.julia/v0.4/Laplacians/src/pcg.jl:57


### chimera
Builds the kth chimeric graph on n vertices. It does this by resetting the random number generator seed. It should captute the state of the generator before that and then return it, but it does not yet.

Builds a chimeric graph on n vertices. The components come from pureRandomGraph, connected by joinGraphs, productGraph and generalizedNecklace


```julia
chimera(n::Integer)
chimera(n::Integer, k::Integer)
```

 at /Users/spielman/.julia/v0.4/Laplacians/src/graphGenerators.jl:515


### compConductance
Returns the quality of the cut for a given graph and a given cut set s.   the result will be |outgoing edges| / min(|vertices in set|, |N - vertices in set|)


```julia
compConductance{Tv,Ti}(G::SparseMatrixCSC{Tv,Ti}, s::Array{Int64,1})
```

 at /Users/spielman/.julia/v0.4/Laplacians/src/graphUtils.jl:163


### compDepth

```julia
compDepth{Tv,Ti}(t::Laplacians.RootedTree{Tv,Ti})
```

 at /Users/spielman/.julia/v0.4/Laplacians/src/treeAlgs.jl:311


### compStretches
Compute the stretched of every edge in `mat` with respect to the tree `tree`. Returns the answer as a sparse matrix with the same nonzero structure as `mat`. Assumes that `mat` is symmetric. `tree` should be the adjacency matrix of a spanning tree.


```julia
compStretches{Tv,Ti}(t::Laplacians.RootedTree{Tv,Ti}, mat::SparseMatrixCSC{Tv,Ti})
compStretches{Tv,Ti}(tree::SparseMatrixCSC{Tv,Ti}, mat::SparseMatrixCSC{Tv,Ti})
```

 at /Users/spielman/.julia/v0.4/Laplacians/src/treeAlgs.jl:393


### completeBinaryTree
The complete binary tree on n vertices


```julia
completeBinaryTree(n::Int64)
```

 at /Users/spielman/.julia/v0.4/Laplacians/src/graphGenerators.jl:139


### completeGraph
The complete graph


```julia
completeGraph(n::Int64)
```

 at /Users/spielman/.julia/v0.4/Laplacians/src/graphGenerators.jl:17


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

 at /Users/spielman/.julia/v0.4/Laplacians/src/graphAlgs.jl:65


### deg

```julia
deg{Tv,Ti}(mat::SparseMatrixCSC{Tv,Ti}, v::Ti)
```

 at /Users/spielman/.julia/v0.4/Laplacians/src/graphUtils.jl:11


### diagmat
Returns the diagonal matrix(as a sparse matrix) of a graph


```julia
diagmat{Tv,Ti}(G::SparseMatrixCSC{Tv,Ti})
```

 at /Users/spielman/.julia/v0.4/Laplacians/src/graphOps.jl:212


### dirEdgeVertexMat
The signed edge-vertex adjacency matrix


```julia
dirEdgeVertexMat(A::SparseMatrixCSC{Tv,Ti<:Integer})
```

 at /Users/spielman/.julia/v0.4/Laplacians/src/toposort.jl:49


### disjoin
Create a disjoint union of graphs a and b,   with no edges between them.


```julia
disjoin(a, b)
```

 at /Users/spielman/.julia/v0.4/Laplacians/src/graphOps.jl:130


### dumb
```
Modify a cluster by passing through all the vertices exactly once and 
adding/removing them based on the value of (Deg_external - Deg_Internal).
```


```julia
dumb{Tv,Ti}(G::SparseMatrixCSC{Tv,Ti}, s::Array{Int64,1})
```

 at /Users/spielman/.julia/v0.4/Laplacians/src/cutHeuristics.jl:106


### edgeVertexMat
The signed edge-vertex adjacency matrix


```julia
edgeVertexMat(mat::SparseMatrixCSC{Tv,Ti<:Integer})
```

 at /Users/spielman/.julia/v0.4/Laplacians/src/graphOps.jl:67


### findEntries
Similar to findnz, but also returns 0 entries that have an edge in the sparse matrix 


```julia
findEntries{Tv,Ti}(G::SparseMatrixCSC{Tv,Ti})
```

 at /Users/spielman/.julia/v0.4/Laplacians/src/graphUtils.jl:141


### flipIndex
For a symmetric matrix, this gives the correspondance between pairs of entries in an ijv. So, ai[ind] = aj[flip[ind]].  For example, 

```
(ai,aj,av) = findnz(a);
fl = flipIndex(a)
ind = 10
@show backind = fl[10]
@show [ai[ind], aj[ind], av[ind]]
@show [ai[backind], aj[backind], av[backind]];

backind = fl[10] = 4
[ai[ind],aj[ind],av[ind]] = [2.0,4.0,0.7]
[ai[backind],aj[backind],av[backind]] = [4.0,2.0,0.7]
```


```julia
flipIndex{Tval,Tind}(a::SparseMatrixCSC{Tval,Tind})
```

 at /Users/spielman/.julia/v0.4/Laplacians/src/graphUtils.jl:51


### floatGraph
Convert the nonzero entries in a graph to Float64


```julia
floatGraph(a::SparseMatrixCSC{Tv,Ti<:Integer})
```

 at /Users/spielman/.julia/v0.4/Laplacians/src/graphOps.jl:6


### generalizedNecklace
Constructs a generalized necklace graph starting with two graphs A and H. The resulting new graph will be constructed by expanding each vertex in H to an instance of A. k random edges will be generated between components. Thus, the resulting graph may have weighted edges.


```julia
generalizedNecklace{Tv,Ti}(A::SparseMatrixCSC{Tv,Ti}, H::SparseMatrixCSC{Tv,Ti<:Integer}, k::Int64)
```

 at /Users/spielman/.julia/v0.4/Laplacians/src/graphOps.jl:232


### generalizedRing
A generalization of a ring graph. The vertices are integers modulo n. Two are connected if their difference is in gens. For example, 

```
generalizedRing(17, [1 5])
```


```julia
generalizedRing(n::Int64, gens)
```

 at /Users/spielman/.julia/v0.4/Laplacians/src/graphGenerators.jl:38


### getObound
Computes the number of edges leaving s 


```julia
getObound{Tv,Ti}(G::SparseMatrixCSC{Tv,Ti}, s::Array{Int64,1})
```

 at /Users/spielman/.julia/v0.4/Laplacians/src/graphUtils.jl:194


### getVolume
Computes the volume of subset s in an unweighted graph G 


```julia
getVolume{Tv,Ti}(G::SparseMatrixCSC{Tv,Ti}, s::Array{Int64,1})
```

 at /Users/spielman/.julia/v0.4/Laplacians/src/graphUtils.jl:176


### grid2
An n-by-m grid graph.  iostropy is the weighting on edges in one direction.


```julia
grid2(n::Int64)
grid2(n::Int64, m::Int64)
```

 at /Users/spielman/.julia/v0.4/Laplacians/src/graphGenerators.jl:159


### grid2coords
Coordinates for plotting the vertices of the n-by-m grid graph


```julia
grid2coords(n::Int64, m::Int64)
grid2coords(n)
```

 at /Users/spielman/.julia/v0.4/Laplacians/src/graphGenerators.jl:163


### grownGraph
Create a graph on n vertices. For each vertex, give it k edges to randomly chosen prior vertices. This is a variety of a preferential attachment graph.    


```julia
grownGraph(n::Int64, k::Int64)
```

 at /Users/spielman/.julia/v0.4/Laplacians/src/graphGenerators.jl:202


### grownGraphD
Like a grownGraph, but it forces the edges to all be distinct. It starts out with a k+1 clique on the first k vertices


```julia
grownGraphD(n::Int64, k::Int64)
```

 at /Users/spielman/.julia/v0.4/Laplacians/src/graphGenerators.jl:234


### hyperCube
The d dimensional hypercube.  Has 2^d vertices


```julia
hyperCube(d::Int64)
```

 at /Users/spielman/.julia/v0.4/Laplacians/src/graphGenerators.jl:125


### isConnected
Returns true if graph is connected.  Calls components.


```julia
isConnected(mat::SparseMatrixCSC{Tv,Ti<:Integer})
```

 at /Users/spielman/.julia/v0.4/Laplacians/src/graphAlgs.jl:113


### joinGraphs
Create a disjoint union of graphs a and b,  and then put k random edges between them


```julia
joinGraphs{Tval,Tind}(a::SparseMatrixCSC{Tval,Tind}, b::SparseMatrixCSC{Tval,Tind}, k::Integer)
```

 at /Users/spielman/.julia/v0.4/Laplacians/src/graphOps.jl:111


### kruskal
Uses Kruskal's algorithm to compute a minimum (or maximum) spanning tree. Set kind=:max if you want the max spanning tree. It returns it a a graph


```julia
kruskal{Tv,Ti}(mat::SparseMatrixCSC{Tv,Ti})
```

 at /Users/spielman/.julia/v0.4/Laplacians/src/graphAlgs.jl:407


### lap
Create a Laplacian matrix from an adjacency matrix. We might want to do this differently, say by enforcing symmetry


```julia
lap(a)
```

 at /Users/spielman/.julia/v0.4/Laplacians/src/graphOps.jl:12


### lapChol


### lapWrapSolver
Takes a solver for solving nonsingular sdd systems, and returns a solver for solving Laplacian systems. The optional args tol and maxits are not necessarily taken by all solvers.  But, if they are, one can pass them here


```julia
lapWrapSolver(solver)
lapWrapSolver(solver, la::AbstractArray{T,N})
lapWrapSolver(solver, la::AbstractArray{T,N}, b)
```

 at /Users/spielman/.julia/v0.4/Laplacians/src/solvers.jl:178


### localImprove
localImprove{Tv,Ti}(G::SparseMatrixCSC{Tv,Ti}, A::Array{Int64,1}; epsSigma=-1.0, err=1e-10, maxSize = max(G.n, G.m)

The LocalImprove function, from the Orrechia-Zhu paper. Given a graph and an initial set, finds a set of smaller conductance based on the starting set using a localized version of max-flow.

Small discussion: When adding in the neighbors of the initial component, if the resulting  conductance is worse than the initial one,  the algorithm will add more and more vertices until hitting a better conductance. However, if we fix a certain  maximum size for our component,  it might be the case that this new conductance will always be worse than what we had initially. Thus, if we run the algorithm with a small maxSize,  our initial conductance might be the best solution we can raech.

  * G is the given graph, A is the initial set 
  * epsSigma is a measure of the quality of the returning set (the smaller the better). It's defaulted to volume(A) / volume(VA)
  * err is the numerical error considered throughout the algorithm. It's defaulted to 1e-10
  * maxSize is the maximum allowed size for the flow graph at any iteration of the algorithm. It's defaulted to |V|


```julia
localImprove{Tv,Ti}(G::SparseMatrixCSC{Tv,Ti}, A::Array{Int64,1})
```

 at /Users/spielman/.julia/v0.4/Laplacians/src/localClustering.jl:22


### mapweight
Create a new graph that is the same as the original, but with f applied to each nonzero entry of a. For example, to make the weight of every edge uniform in [0,1], we could write

```julia
b = mapweight(a, x->rand(1)[1])
```


```julia
mapweight{Tval,Tind}(a::SparseMatrixCSC{Tval,Tind}, f)
```

 at /Users/spielman/.julia/v0.4/Laplacians/src/graphOps.jl:40


### matToTree

```julia
matToTree{Tv,Ti}(mat::SparseMatrixCSC{Tv,Ti})
matToTree{Tv,Ti}(mat::SparseMatrixCSC{Tv,Ti}, root::Ti)
```

 at /Users/spielman/.julia/v0.4/Laplacians/src/treeAlgs.jl:32


### matToTreeDepth

```julia
matToTreeDepth{Tv,Ti}(mat::SparseMatrixCSC{Tv,Ti})
matToTreeDepth{Tv,Ti}(mat::SparseMatrixCSC{Tv,Ti}, root::Ti)
```

 at /Users/spielman/.julia/v0.4/Laplacians/src/treeAlgs.jl:98


### maxflow
implementation of Dinic's algorithm. computes the maximum flow and min-cut in G between s and t    we consider the adjacency matrix to be the capacity matrix 


```julia
maxflow{Tv,Ti}(G::SparseMatrixCSC{Tv,Ti}, s::Int64, t::Int64)
```

 at /Users/spielman/.julia/v0.4/Laplacians/src/flow.jl:7


### nbri

```julia
nbri{Tv,Ti}(mat::SparseMatrixCSC{Tv,Ti}, v::Ti, i::Ti)
```

 at /Users/spielman/.julia/v0.4/Laplacians/src/graphUtils.jl:12


### nbrs

```julia
nbrs{Tv,Ti}(mat::SparseMatrixCSC{Tv,Ti}, v::Ti)
```

 at /Users/spielman/.julia/v0.4/Laplacians/src/graphUtils.jl:14


### pathFromParents

```julia
pathFromParents(parents, y)
```

 at /Users/spielman/.julia/v0.4/Laplacians/src/graphAlgs.jl:215


### pathGraph
The path graph on n vertices


```julia
pathGraph(n::Int64)
```

 at /Users/spielman/.julia/v0.4/Laplacians/src/graphGenerators.jl:8


### pcg
`pcg(mat, b, pre; tol, maxits, maxtime, verbose)` solves a symmetric linear system using preconditioner `pre`. `pre` can be a function or a matrix.  If a matrix, a function to solve it is created with cholFact. `tol` is set to 1e-6 by default, `maxits` defaults to Inf `maxtime` defaults to Inf.  It measures seconds. `verbose` defaults to false


```julia
pcg(mat, b, pre::AbstractArray{T,N})
pcg(mat, b::Array{Float64,1}, pre::Function)
pcg(mat, b::Array{Float32,1}, pre::Function)
pcg(mat, b, pre::Function)
```

 at /Users/spielman/.julia/v0.4/Laplacians/src/pcg.jl:63


### pcgLapSolver
Create a solver that uses cg to solve Laplacian systems in mat. Specialized for the case when pre is a Laplacian matrix.  Fix the default parameters of the solver as given


```julia
pcgLapSolver(mat, pre)
```

 at /Users/spielman/.julia/v0.4/Laplacians/src/pcg.jl:91


### pcgSolver
`pcgSolver(mat, pre; tol, maxits, maxtime, verbose)` creates a solver for a PSD system using preconditioner `pre`. The parameters are as described in pcg.


```julia
pcgSolver(mat, pre)
```

 at /Users/spielman/.julia/v0.4/Laplacians/src/pcg.jl:85


### plotGraph
Plots graph gr with coordinates (x,y)


```julia
plotGraph(gr, x, y)
plotGraph(gr, x, y, color)
```

 at /Users/spielman/.julia/v0.4/Laplacians/src/graphOps.jl:136


### prefAttach
A preferential attachment graph in which each vertex has k edges to those that come before.  These are chosen with probability p to be from a random vertex, and with probability 1-p to come from the endpoint of a random edge. It begins with a k-clique on the first k+1 vertices.


```julia
prefAttach(n::Int64, k::Int64, p::Float64)
```

 at /Users/spielman/.julia/v0.4/Laplacians/src/graphGenerators.jl:258


### prim
`prim(mat::SparseMatrixCSC; rev=false)` Compute a minimum spanning tree of the matrix `mat`.   If rev is true, computes a maximum spanning tree.


```julia
prim(mat::SparseMatrixCSC{Tv,Ti<:Integer})
```

 at /Users/spielman/.julia/v0.4/Laplacians/src/graphAlgs.jl:440


### prn
prn{Tv, Ti}(G::SparseMatrixCSC{Tv,Ti}, s::Array{Int64,1}, phi::Float64, b::Int64)

The PageRank-Nibble cutting algorithm from the Anderson/Chung/Lang paper

s is a set of starting vertices, phi is a constant in (0, 1], and b is an integer in [1, [log m]]

phi is a bound on the quality of the conductance of the cut - the smaller the phi, the higher the quality.  b is used to handle precision throughout the algorithm - the higher the b, the greater the precision.


```julia
prn{Tv,Ti}(G::SparseMatrixCSC{Tv,Ti}, s::Array{Int64,1}, phi::Float64, b::Int64)
```

 at /Users/spielman/.julia/v0.4/Laplacians/src/localClustering.jl:374


### productGraph
The Cartesian product of two graphs.  When applied to two paths, it gives a grid.


```julia
productGraph(a0::SparseMatrixCSC{Tv,Ti<:Integer}, a1::SparseMatrixCSC{Tv,Ti<:Integer})
```

 at /Users/spielman/.julia/v0.4/Laplacians/src/graphOps.jl:58


### pureRandomGraph
Generate a random graph with n vertices from one of our natural distributions


```julia
pureRandomGraph(n::Integer)
```

 at /Users/spielman/.julia/v0.4/Laplacians/src/graphGenerators.jl:379


### randGenRing
A random generalized ring graph of degree k. Gens always contains 1, and the other k-1 edge types are chosen from an exponential distribution


```julia
randGenRing(n::Int64, k::Integer)
```

 at /Users/spielman/.julia/v0.4/Laplacians/src/graphGenerators.jl:62


### randMatching
A random matching on n vertices


```julia
randMatching(n::Int64)
```

 at /Users/spielman/.julia/v0.4/Laplacians/src/graphGenerators.jl:174


### randRegular
A sum of k random matchings on n vertices


```julia
randRegular(n::Int64, k::Int64)
```

 at /Users/spielman/.julia/v0.4/Laplacians/src/graphGenerators.jl:187


### randWeight
Applies one of a number of random weighting schemes to the edges of the graph


```julia
randWeight(a)
```

 at /Users/spielman/.julia/v0.4/Laplacians/src/graphGenerators.jl:535


### randishKruskal
A heuristic for computing low-stretch spanning trees.  Where Kruskal's MST algorithm adds edges in order of weight, this algorithm adds them at random with probability proportional to their weight.


```julia
randishKruskal{Tv,Ti}(mat::SparseMatrixCSC{Tv,Ti})
```

 at /Users/spielman/.julia/v0.4/Laplacians/src/randTrees.jl:11


### randishPrim
A heuristic for computing low-stretch spanning trees.  Where Prim's MST algorithm grows a cluster by always adding the edge on the boundary of maximum weight, this algorithm adds a boundary edge with probability proportional to its weight.


```julia
randishPrim{Tval,Tind}(mat::SparseMatrixCSC{Tval,Tind})
```

 at /Users/spielman/.julia/v0.4/Laplacians/src/randTrees.jl:40


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

 at /Users/spielman/.julia/v0.4/Laplacians/src/IO.jl:4


### readIJV
To read a simple edge list, each line being an (i, j, v) pair. The parens should not be there in the format, just commas separating. To generate this format in Matlab, you just need to be careful to write the vertex indices with sufficient precision.  For example, you can do this

```
>> [ai,aj,av] = find(triu(a));
>> dlmwrite('graph.txt',[ai,aj,av],'precision',9);
```


```julia
readIJV(filename::AbstractString)
```

 at /Users/spielman/.julia/v0.4/Laplacians/src/IO.jl:25


### refineCut
```
Modify a cluster by adding or removing vertices by picking at each step 
the vertex that has the maximum value of (Deg_external - Deg_Internal).
Each vertex can be added in/removed only once.
```


```julia
refineCut{Tv,Ti}(G::SparseMatrixCSC{Tv,Ti}, s::Array{Int64,1})
```

 at /Users/spielman/.julia/v0.4/Laplacians/src/cutHeuristics.jl:9


### ringGraph
The simple ring on n vertices


```julia
ringGraph(n::Int64)
```

 at /Users/spielman/.julia/v0.4/Laplacians/src/graphGenerators.jl:23


### semiWtedChimera
A Chimera graph with some weights.  The weights just appear when graphs are combined. For more interesting weights, use `wtedChimera`


```julia
semiWtedChimera(n::Integer)
semiWtedChimera(n::Integer, k::Integer)
```

 at /Users/spielman/.julia/v0.4/Laplacians/src/graphGenerators.jl:442


### setValue
Sets the value of a certain edge in a sparse graph; value can be 0 without the edges dissapearing 


```julia
setValue{Tv,Ti}(mat::SparseMatrixCSC{Tv,Ti}, v::Ti, i::Ti, a::Tv)
```

 at /Users/spielman/.julia/v0.4/Laplacians/src/graphUtils.jl:29


### shortIntGraph
Convert the indices in a graph to 32-bit ints.  This takes less storage, but does not speed up much


```julia
shortIntGraph(a::SparseMatrixCSC{Tv,Ti<:Integer})
```

 at /Users/spielman/.julia/v0.4/Laplacians/src/graphOps.jl:3


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

 at /Users/spielman/.julia/v0.4/Laplacians/src/graphAlgs.jl:239


### shortestPaths
Computes the lenghts of shortest paths from `start`. Returns both a vector of the lenghts, and the parent array in the shortest path tree.

This algorithm treats edge weights as reciprocals of distances. DOC BETTER


```julia
shortestPaths{Tv,Ti}(mat::SparseMatrixCSC{Tv,Ti}, start::Ti)
shortestPaths{Tv,Ti}(mat::SparseMatrixCSC{Tv,Ti})
```

 at /Users/spielman/.julia/v0.4/Laplacians/src/graphAlgs.jl:175


### spectralCoords
Computes the spectral coordinates of a graph


```julia
spectralCoords(a)
```

 at /Users/spielman/.julia/v0.4/Laplacians/src/graphOps.jl:181


### spectralDrawing
Computes spectral coordinates, and then uses plotGraph to draw


```julia
spectralDrawing(a)
```

 at /Users/spielman/.julia/v0.4/Laplacians/src/graphOps.jl:173


### subsampleEdges
Create a new graph from the old, but keeping edge edge with probability `p`


```julia
subsampleEdges(a::SparseMatrixCSC{Float64,Int64}, p::Float64)
```

 at /Users/spielman/.julia/v0.4/Laplacians/src/graphOps.jl:76


### tarjanStretch

```julia
tarjanStretch{Tv,Ti}(t::Laplacians.RootedTree{Tv,Ti}, mat::SparseMatrixCSC{Tv,Ti}, depth::Array{Tv,1})
```

 at /Users/spielman/.julia/v0.4/Laplacians/src/treeAlgs.jl:334


### toUnitVector
Creates a unit vector of length n from a given set of integers, with weights based on the number of occurences


```julia
toUnitVector(a::Array{Int64,1}, n)
```

 at /Users/spielman/.julia/v0.4/Laplacians/src/graphOps.jl:190


### toposort

```julia
toposort{Tv,Ti}(mat::SparseMatrixCSC{Tv,Ti})
```

 at /Users/spielman/.julia/v0.4/Laplacians/src/toposort.jl:13


### twoLift
Creats a 2-lift of a.  `flip` is a boolean indicating which edges cross


```julia
twoLift(a)
twoLift(a, flip::AbstractArray{Bool,1})
twoLift(a, k::Integer)
```

 at /Users/spielman/.julia/v0.4/Laplacians/src/graphOps.jl:99


### uniformWeight
Put a uniform [0,1] weight on every edge.  This is an example of how to use mapweight.


```julia
uniformWeight{Tval,Tind}(a::SparseMatrixCSC{Tval,Tind})
```

 at /Users/spielman/.julia/v0.4/Laplacians/src/graphOps.jl:49


### uniformWeight!
Set the weight of every edge to 1


```julia
uniformWeight!(mat::SparseMatrixCSC{Tv,Ti<:Integer})
```

 at /Users/spielman/.julia/v0.4/Laplacians/src/graphOps.jl:53


### unweight
Create a new graph in that is the same as the original, but with all edge weights 1


```julia
unweight{Tval,Tind}(ain::SparseMatrixCSC{Tval,Tind})
```

 at /Users/spielman/.julia/v0.4/Laplacians/src/graphOps.jl:16


### unweight!
Change the weight of every edge in a to 1


```julia
unweight!{Tval,Tind}(a::SparseMatrixCSC{Tval,Tind})
```

 at /Users/spielman/.julia/v0.4/Laplacians/src/graphOps.jl:26


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

 at /Users/spielman/.julia/v0.4/Laplacians/src/graphAlgs.jl:136


### wdeg
Finds the weighted degree of a vertex in the graph 


```julia
wdeg{Tv,Ti}(mat::SparseMatrixCSC{Tv,Ti}, v::Ti)
```

 at /Users/spielman/.julia/v0.4/Laplacians/src/graphUtils.jl:19


### weighti

```julia
weighti{Tv,Ti}(mat::SparseMatrixCSC{Tv,Ti}, v::Ti, i::Ti)
```

 at /Users/spielman/.julia/v0.4/Laplacians/src/graphUtils.jl:13


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

 at /Users/spielman/.julia/v0.4/Laplacians/src/IO.jl:52


### wtedChimera
Builds the kth wted chimeric graph on n vertices. It does this by resetting the random number generator seed. It should captute the state of the generator before that and then return it, but it does not yet.

Generate a chimera, and then apply a random weighting scheme


```julia
wtedChimera(n::Integer)
wtedChimera(n::Integer, k::Integer)
```

 at /Users/spielman/.julia/v0.4/Laplacians/src/graphGenerators.jl:607


