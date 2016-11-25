
<a id='graph-Generators-1'></a>

# graph Generators


*Laplacians.jl* implements generators for many standard graphs. It the `chimera` and `wtedChimera` generators combine these standard graphs by


  * joining disjoint copies with edges,
  * forming Kronecker products, and
  * taking `generalizedNecklace` products.


<a id='Base.Random.randperm-Tuple{AbstractArray{T,2}}' href='#Base.Random.randperm-Tuple{AbstractArray{T,2}}'>#</a>
**`Base.Random.randperm`** &mdash; *Method*.



Randomly permutes the vertex indices


<a target='_blank' href='https://github.com/danspielman/Laplacians.jl/tree/ce307408e6b2290f943763f162deee084e2b9097/doc/../src/graphGenerators.jl#L313' class='documenter-source'>source</a><br>

<a id='Laplacians.ErdosRenyi-Tuple{Integer,Integer}' href='#Laplacians.ErdosRenyi-Tuple{Integer,Integer}'>#</a>
**`Laplacians.ErdosRenyi`** &mdash; *Method*.



Generate a random graph on n vertices with m edges. The actual number of edges will probably be smaller, as we sample with replacement


<a target='_blank' href='https://github.com/danspielman/Laplacians.jl/tree/ce307408e6b2290f943763f162deee084e2b9097/doc/../src/graphGenerators.jl#L322-L325' class='documenter-source'>source</a><br>

<a id='Laplacians.ErdosRenyiCluster-Tuple{Integer,Integer}' href='#Laplacians.ErdosRenyiCluster-Tuple{Integer,Integer}'>#</a>
**`Laplacians.ErdosRenyiCluster`** &mdash; *Method*.



Generate an ER graph with average degree k, and then return the largest component. Will probably have fewer than n vertices. If you want to add a tree to bring it back to n, try ErdosRenyiClusterFix.


<a target='_blank' href='https://github.com/danspielman/Laplacians.jl/tree/ce307408e6b2290f943763f162deee084e2b9097/doc/../src/graphGenerators.jl#L335-L340' class='documenter-source'>source</a><br>

<a id='Laplacians.ErdosRenyiClusterFix-Tuple{Integer,Integer}' href='#Laplacians.ErdosRenyiClusterFix-Tuple{Integer,Integer}'>#</a>
**`Laplacians.ErdosRenyiClusterFix`** &mdash; *Method*.



Like an Erdos-Renyi cluster, but add back a tree so it has n vertices


<a target='_blank' href='https://github.com/danspielman/Laplacians.jl/tree/ce307408e6b2290f943763f162deee084e2b9097/doc/../src/graphGenerators.jl#L351-L353' class='documenter-source'>source</a><br>

<a id='Laplacians.chimera-Tuple{Integer,Integer}' href='#Laplacians.chimera-Tuple{Integer,Integer}'>#</a>
**`Laplacians.chimera`** &mdash; *Method*.



Builds the kth chimeric graph on n vertices. It does this by resetting the random number generator seed. It should captute the state of the generator before that and then return it, but it does not yet.


<a target='_blank' href='https://github.com/danspielman/Laplacians.jl/tree/ce307408e6b2290f943763f162deee084e2b9097/doc/../src/graphGenerators.jl#L510-L514' class='documenter-source'>source</a><br>

<a id='Laplacians.chimera-Tuple{Integer}' href='#Laplacians.chimera-Tuple{Integer}'>#</a>
**`Laplacians.chimera`** &mdash; *Method*.



Builds a chimeric graph on n vertices. The components come from pureRandomGraph, connected by joinGraphs, productGraph and generalizedNecklace


<a target='_blank' href='https://github.com/danspielman/Laplacians.jl/tree/ce307408e6b2290f943763f162deee084e2b9097/doc/../src/graphGenerators.jl#L497-L500' class='documenter-source'>source</a><br>

<a id='Laplacians.completeBinaryTree-Tuple{Int64}' href='#Laplacians.completeBinaryTree-Tuple{Int64}'>#</a>
**`Laplacians.completeBinaryTree`** &mdash; *Method*.



The complete binary tree on n vertices


<a target='_blank' href='https://github.com/danspielman/Laplacians.jl/tree/ce307408e6b2290f943763f162deee084e2b9097/doc/../src/graphGenerators.jl#L85' class='documenter-source'>source</a><br>

<a id='Laplacians.completeGraph-Tuple{Int64}' href='#Laplacians.completeGraph-Tuple{Int64}'>#</a>
**`Laplacians.completeGraph`** &mdash; *Method*.



The complete graph


<a target='_blank' href='https://github.com/danspielman/Laplacians.jl/tree/ce307408e6b2290f943763f162deee084e2b9097/doc/../src/graphGenerators.jl#L15' class='documenter-source'>source</a><br>

<a id='Laplacians.generalizedRing-Tuple{Int64,Any}' href='#Laplacians.generalizedRing-Tuple{Int64,Any}'>#</a>
**`Laplacians.generalizedRing`** &mdash; *Method*.



A generalization of a ring graph. The vertices are integers modulo n. Two are connected if their difference is in gens. For example, 

```
generalizedRing(17, [1 5])
```


<a target='_blank' href='https://github.com/danspielman/Laplacians.jl/tree/ce307408e6b2290f943763f162deee084e2b9097/doc/../src/graphGenerators.jl#L28-L37' class='documenter-source'>source</a><br>

<a id='Laplacians.grid2-Tuple{Int64,Int64}' href='#Laplacians.grid2-Tuple{Int64,Int64}'>#</a>
**`Laplacians.grid2`** &mdash; *Method*.



An n-by-m grid graph.  iostropy is the weighting on edges in one direction.


<a target='_blank' href='https://github.com/danspielman/Laplacians.jl/tree/ce307408e6b2290f943763f162deee084e2b9097/doc/../src/graphGenerators.jl#L132' class='documenter-source'>source</a><br>

<a id='Laplacians.grid2coords-Tuple{Int64,Int64}' href='#Laplacians.grid2coords-Tuple{Int64,Int64}'>#</a>
**`Laplacians.grid2coords`** &mdash; *Method*.



Coordinates for plotting the vertices of the n-by-m grid graph


<a target='_blank' href='https://github.com/danspielman/Laplacians.jl/tree/ce307408e6b2290f943763f162deee084e2b9097/doc/../src/graphGenerators.jl#L150' class='documenter-source'>source</a><br>

<a id='Laplacians.grid3-Tuple{Ti,Ti,Ti}' href='#Laplacians.grid3-Tuple{Ti,Ti,Ti}'>#</a>
**`Laplacians.grid3`** &mdash; *Method*.



An n1-by-n2-by-n3 grid graph.


<a target='_blank' href='https://github.com/danspielman/Laplacians.jl/tree/ce307408e6b2290f943763f162deee084e2b9097/doc/../src/graphGenerators.jl#L142' class='documenter-source'>source</a><br>

<a id='Laplacians.grownGraph-Tuple{Int64,Int64}' href='#Laplacians.grownGraph-Tuple{Int64,Int64}'>#</a>
**`Laplacians.grownGraph`** &mdash; *Method*.



Create a graph on n vertices. For each vertex, give it k edges to randomly chosen prior vertices. This is a variety of a preferential attachment graph.    


<a target='_blank' href='https://github.com/danspielman/Laplacians.jl/tree/ce307408e6b2290f943763f162deee084e2b9097/doc/../src/graphGenerators.jl#L185-L190' class='documenter-source'>source</a><br>

<a id='Laplacians.grownGraphD-Tuple{Int64,Int64}' href='#Laplacians.grownGraphD-Tuple{Int64,Int64}'>#</a>
**`Laplacians.grownGraphD`** &mdash; *Method*.



Like a grownGraph, but it forces the edges to all be distinct. It starts out with a k+1 clique on the first k vertices


<a target='_blank' href='https://github.com/danspielman/Laplacians.jl/tree/ce307408e6b2290f943763f162deee084e2b9097/doc/../src/graphGenerators.jl#L220-L222' class='documenter-source'>source</a><br>

<a id='Laplacians.hyperCube-Tuple{Int64}' href='#Laplacians.hyperCube-Tuple{Int64}'>#</a>
**`Laplacians.hyperCube`** &mdash; *Method*.



The d dimensional hypercube.  Has 2^d vertices


<a target='_blank' href='https://github.com/danspielman/Laplacians.jl/tree/ce307408e6b2290f943763f162deee084e2b9097/doc/../src/graphGenerators.jl#L72' class='documenter-source'>source</a><br>

<a id='Laplacians.pathGraph-Tuple{Int64}' href='#Laplacians.pathGraph-Tuple{Int64}'>#</a>
**`Laplacians.pathGraph`** &mdash; *Method*.



The path graph on n vertices


<a target='_blank' href='https://github.com/danspielman/Laplacians.jl/tree/ce307408e6b2290f943763f162deee084e2b9097/doc/../src/graphGenerators.jl#L6' class='documenter-source'>source</a><br>

<a id='Laplacians.prefAttach-Tuple{Int64,Int64,Float64}' href='#Laplacians.prefAttach-Tuple{Int64,Int64,Float64}'>#</a>
**`Laplacians.prefAttach`** &mdash; *Method*.



A preferential attachment graph in which each vertex has k edges to those that come before.  These are chosen with probability p to be from a random vertex, and with probability 1-p to come from the endpoint of a random edge. It begins with a k-clique on the first k+1 vertices.


<a target='_blank' href='https://github.com/danspielman/Laplacians.jl/tree/ce307408e6b2290f943763f162deee084e2b9097/doc/../src/graphGenerators.jl#L242-L246' class='documenter-source'>source</a><br>

<a id='Laplacians.pureRandomGraph-Tuple{Integer}' href='#Laplacians.pureRandomGraph-Tuple{Integer}'>#</a>
**`Laplacians.pureRandomGraph`** &mdash; *Method*.



Generate a random graph with n vertices from one of our natural distributions


<a target='_blank' href='https://github.com/danspielman/Laplacians.jl/tree/ce307408e6b2290f943763f162deee084e2b9097/doc/../src/graphGenerators.jl#L366' class='documenter-source'>source</a><br>

<a id='Laplacians.randGenRing-Tuple{Int64,Integer}' href='#Laplacians.randGenRing-Tuple{Int64,Integer}'>#</a>
**`Laplacians.randGenRing`** &mdash; *Method*.



A random generalized ring graph of degree k. Gens always contains 1, and the other k-1 edge types are chosen from an exponential distribution


<a target='_blank' href='https://github.com/danspielman/Laplacians.jl/tree/ce307408e6b2290f943763f162deee084e2b9097/doc/../src/graphGenerators.jl#L57-L60' class='documenter-source'>source</a><br>

<a id='Laplacians.randMatching-Tuple{Int64}' href='#Laplacians.randMatching-Tuple{Int64}'>#</a>
**`Laplacians.randMatching`** &mdash; *Method*.



A random matching on n vertices


<a target='_blank' href='https://github.com/danspielman/Laplacians.jl/tree/ce307408e6b2290f943763f162deee084e2b9097/doc/../src/graphGenerators.jl#L160' class='documenter-source'>source</a><br>

<a id='Laplacians.randRegular-Tuple{Int64,Int64}' href='#Laplacians.randRegular-Tuple{Int64,Int64}'>#</a>
**`Laplacians.randRegular`** &mdash; *Method*.



A sum of k random matchings on n vertices


<a target='_blank' href='https://github.com/danspielman/Laplacians.jl/tree/ce307408e6b2290f943763f162deee084e2b9097/doc/../src/graphGenerators.jl#L174' class='documenter-source'>source</a><br>

<a id='Laplacians.randWeight-Tuple{Any}' href='#Laplacians.randWeight-Tuple{Any}'>#</a>
**`Laplacians.randWeight`** &mdash; *Method*.



Applies one of a number of random weighting schemes to the edges of the graph


<a target='_blank' href='https://github.com/danspielman/Laplacians.jl/tree/ce307408e6b2290f943763f162deee084e2b9097/doc/../src/graphGenerators.jl#L520' class='documenter-source'>source</a><br>

<a id='Laplacians.ringGraph-Tuple{Int64}' href='#Laplacians.ringGraph-Tuple{Int64}'>#</a>
**`Laplacians.ringGraph`** &mdash; *Method*.



The simple ring on n vertices


<a target='_blank' href='https://github.com/danspielman/Laplacians.jl/tree/ce307408e6b2290f943763f162deee084e2b9097/doc/../src/graphGenerators.jl#L21' class='documenter-source'>source</a><br>

<a id='Laplacians.semiWtedChimera-Tuple{Integer}' href='#Laplacians.semiWtedChimera-Tuple{Integer}'>#</a>
**`Laplacians.semiWtedChimera`** &mdash; *Method*.



A Chimera graph with some weights.  The weights just appear when graphs are combined. For more interesting weights, use `wtedChimera`


<a target='_blank' href='https://github.com/danspielman/Laplacians.jl/tree/ce307408e6b2290f943763f162deee084e2b9097/doc/../src/graphGenerators.jl#L428-L430' class='documenter-source'>source</a><br>

<a id='Laplacians.wGrid2-Tuple{Int64}' href='#Laplacians.wGrid2-Tuple{Int64}'>#</a>
**`Laplacians.wGrid2`** &mdash; *Method*.



An n by n grid with random weights. User can specify the weighting scheme. 


<a target='_blank' href='https://github.com/danspielman/Laplacians.jl/tree/ce307408e6b2290f943763f162deee084e2b9097/doc/../src/graphGenerators.jl#L100' class='documenter-source'>source</a><br>

<a id='Laplacians.wGrid3-Tuple{Int64}' href='#Laplacians.wGrid3-Tuple{Int64}'>#</a>
**`Laplacians.wGrid3`** &mdash; *Method*.



An n^3 grid with random weights. User can specify the weighting scheme. 


<a target='_blank' href='https://github.com/danspielman/Laplacians.jl/tree/ce307408e6b2290f943763f162deee084e2b9097/doc/../src/graphGenerators.jl#L114' class='documenter-source'>source</a><br>

<a id='Laplacians.wtedChimera-Tuple{Integer,Integer}' href='#Laplacians.wtedChimera-Tuple{Integer,Integer}'>#</a>
**`Laplacians.wtedChimera`** &mdash; *Method*.



Builds the kth wted chimeric graph on n vertices. It does this by resetting the random number generator seed. It should captute the state of the generator before that and then return it, but it does not yet.


<a target='_blank' href='https://github.com/danspielman/Laplacians.jl/tree/ce307408e6b2290f943763f162deee084e2b9097/doc/../src/graphGenerators.jl#L576-L580' class='documenter-source'>source</a><br>

<a id='Laplacians.wtedChimera-Tuple{Integer}' href='#Laplacians.wtedChimera-Tuple{Integer}'>#</a>
**`Laplacians.wtedChimera`** &mdash; *Method*.



Generate a chimera, and then apply a random weighting scheme


<a target='_blank' href='https://github.com/danspielman/Laplacians.jl/tree/ce307408e6b2290f943763f162deee084e2b9097/doc/../src/graphGenerators.jl#L593' class='documenter-source'>source</a><br>

