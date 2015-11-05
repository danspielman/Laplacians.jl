# graphGenerators
### pathGraph
The path graph on n vertices


```julia
pathGraph(n::Int64)
```

graphGenerators.jl:45



### completeGraph
The complete graph


```julia
completeGraph(n::Int64)
```

graphGenerators.jl:54



### ringGraph
The simple ring on n vertices


```julia
ringGraph(n::Int64)
```

graphGenerators.jl:60



### generalizedRing
A generalization of a ring graph. The vertices are integers modulo n. Two are connected if their difference is in gens. For example, 

```
generalizedRing(17, [1 5])
```


```julia
generalizedRing(n::Int64, gens)
```

graphGenerators.jl:75



### randGenRing
A random generalized ring graph of degree k. Gens always contains 1, and the other k-1 edge types are chosen from an exponential distribution


```julia
randGenRing(n::Int64, k::Integer)
```

graphGenerators.jl:98



### hyperCube
The d dimensional hypercube.  Has 2^d vertices


```julia
hyperCube(d::Int64)
```

graphGenerators.jl:105



### completeBinaryTree
The complete binary tree on n vertices


```julia
completeBinaryTree(n::Int64)
```

graphGenerators.jl:119



### grid2
An n-by-m grid graph.  iostropy is the weighting on edges in one direction.


```julia
grid2(n::Int64)
grid2(n::Int64, m::Int64)
```

graphGenerators.jl:139



### grid2coords
Coordinates for plotting the vertices of the n-by-m grid graph


```julia
grid2coords(n::Int64, m::Int64)
grid2coords(n)
```

graphGenerators.jl:143



### randMatching
A random matching on n vertices


```julia
randMatching(n::Int64)
```

graphGenerators.jl:154



### randRegular
A sum of k random matchings on n vertices


```julia
randRegular(n::Int64, k::Int64)
```

graphGenerators.jl:167



### grownGraph
Create a graph on n vertices. For each vertex, give it k edges to randomly chosen prior vertices. This is a variety of a preferential attachment graph.    


```julia
grownGraph(n::Int64, k::Int64)
```

graphGenerators.jl:182



### grownGraphD
Like a grownGraph, but it forces the edges to all be distinct. It starts out with a k+1 clique on the first k vertices


```julia
grownGraphD(n::Int64, k::Int64)
```

graphGenerators.jl:214



### prefAttach
A preferential attachment graph in which each vertex has k edges to those that come before.  These are chosen with probability p to be from a random vertex, and with probability 1-p to come from the endpoint of a random edge. It begins with a k-clique on the first k+1 vertices.


```julia
prefAttach(n::Int64, k::Int64, p::Float64)
```

graphGenerators.jl:238



### ErdosRenyiCluster
Generate an ER graph with average degree k, and then return the largest component. Will probably have fewer than n vertices. If you want to add a tree to bring it back to n, try ErdosRenyiClusterFix.


```julia
ErdosRenyiCluster(n::Integer, k::Integer)
```

graphGenerators.jl:319



### ErdosRenyiClusterFix
Like an Erdos-Renyi cluster, but add back a tree so it has n vertices


```julia
ErdosRenyiClusterFix(n::Integer, k::Integer)
```

graphGenerators.jl:332



### pureRandomGraph
Generate a random graph with n vertices from one of our natural distributions


```julia
pureRandomGraph(n::Integer)
```

graphGenerators.jl:347



### chimera
Builds the kth chimeric graph on n vertices. It does this by resetting the random number generator seed. It should captute the state of the generator before that and then return it, but it does not yet.

Builds a chimeric graph on n vertices. The components come from pureRandomGraph, connected by joinGraphs, productGraph and generalizedNecklace


```julia
chimera(n::Integer)
chimera(n::Integer, k::Integer)
```

graphGenerators.jl:392



### randWeight
Applies one of a number of random weighting schemes to the edges of the graph


```julia
randWeight(a)
```

graphGenerators.jl:463



### wtedChimera
Builds the kth wted chimeric graph on n vertices. It does this by resetting the random number generator seed. It should captute the state of the generator before that and then return it, but it does not yet.

Generate a chimera, and then apply a random weighting scheme


```julia
wtedChimera(n::Integer)
wtedChimera(n::Integer, k::Integer)
```

graphGenerators.jl:528



