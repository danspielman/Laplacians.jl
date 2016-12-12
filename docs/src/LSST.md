# Low Stretch Spanning Trees

We have implemented a variant of the algorithm of Alon, Karp, Peleg and West for computing low stretch spanning trees.  It is called `akpw`.  For unweighted graphs we provide a faster variant called `akpwU`.  If you require faster algorithms, with possibly higher average stretch, you can try the heuristics `randishPrim` or `randishKruskal`.

You can compute the stretch of a graph with respect to a spanning tree with the routine `compStretches`.  It returns a sparse matrix with one entry for each edge in the graph giving its stretch.  For example:

~~~julia
julia> graph = grid2(4)
16x16 sparse matrix with 48 Float64 entries:

julia> tree = akpwU(graph)
16x16 sparse matrix with 30 Float64 entries:

julia> st = compStretches(tree,graph)
16x16 sparse matrix with 48 Float64 entries:
	[2 ,  1]  =  1.0
	[5 ,  1]  =  1.0
	[1 ,  2]  =  1.0
	[3 ,  2]  =  1.0
	[6 ,  2]  =  3.0
	[2 ,  3]  =  1.0
	[4 ,  3]  =  1.0
	[7 ,  3]  =  1.0
	[3 ,  4]  =  1.0
	[8 ,  4]  =  3.0
	[1 ,  5]  =  1.0
	[6 ,  5]  =  5.0
	⋮
	[8 , 12]  =  3.0
	[11, 12]  =  1.0
	[16, 12]  =  1.0
	[9 , 13]  =  1.0
	[14, 13]  =  3.0
	[10, 14]  =  1.0
	[13, 14]  =  3.0
	[15, 14]  =  3.0
	[11, 15]  =  1.0
	[14, 15]  =  3.0
	[16, 15]  =  3.0
	[12, 16]  =  1.0
	[15, 16]  =  3.0
~~~

Here is an example demonstrating the average stretches and times taken by these algorithms on a large graph.

~~~julia

julia> graph = chimera(1000000,1);

julia> @time tree = akpw(graph);
  5.700285 seconds (16.10 M allocations: 1.263 GB, 11.16% gc time)

julia> aveStretch = sum(compStretches(tree,graph))/nnz(graph)
8.793236275779616

julia> @time tree = randishPrim(graph);
  3.736225 seconds (3.21 M allocations: 566.887 MB, 6.40% gc time)

julia> aveStretch = sum(compStretches(tree,graph))/nnz(graph)
10.800094649795756

julia> @time tree = randishKruskal(graph);
  2.515443 seconds (3.21 M allocations: 423.529 MB, 4.35% gc time)

julia> aveStretch = sum(compStretches(tree,graph))/nnz(graph)
37.819948689847564

~~~

Of course, you can get very different results on very different graphs.  But, the ordering of these algorithms respect to time and stretch will usually follow this pattern.
