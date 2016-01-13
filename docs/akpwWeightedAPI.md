# akpwWeighted
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

akpwWeighted.jl:733



### akpw
This is a wrapper for akpw!. It's slower, but won't modify the original graph. See akpw! documentation for more details.


```julia
akpw(origMat::SparseMatrixCSC{Tv,Ti<:Integer})
```

akpwWeighted.jl:818



