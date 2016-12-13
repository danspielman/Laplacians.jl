# Using sparse matrices as graphs

The routines `deg`, `nbri` and `weighti` will let you treat a sparse matrix like a graph.

`deg(graph, u)` is the degree of node u.
`nbri(graph, u, i)` is the ith neighbor of node u.
`weighti(graph, u, i)` is the weight of the edge to the ith neighbor of node u.

Note that we start indexing from 1.

For example, to iterate over the neighbors of node v,
  and play with the attached nodes, you could write code like:

~~~julia
  for i in 1:deg(mat, v)
     nbr = nbri(mat, v, i)
     wt = weighti(mat, v, i)
     foo(v, nbr, wt)
  end
~~~

But, this turns out to be much slower than working with the structure directly, like

~~~julia
  for ind in mat.colptr[v]:(mat.colptr[v+1]-1)
      nbr = mat.rowval[ind]
      wt = mat.nzval[ind]
      foo(v, nbr, wt)
  end
~~~

* [ ] Maybe we can make a macro to replace those functions.  It could be faster and more readable.

### The SparseMatrixCSC data structure

You can explore what is going on with the data structure by looking at some examples.  For example, here is a randomly weighted complete graph on 4 vertices, first displayed as a matrix:

~~~julia
gr = round(10*uniformWeight(completeGraph(4)))

4x4 sparse matrix with 12 Float64 entries:
	[2, 1]  =  3.0
	[3, 1]  =  3.0
	[4, 1]  =  6.0
	[1, 2]  =  3.0
	[3, 2]  =  1.0
	[4, 2]  =  2.0
	[1, 3]  =  3.0
	[2, 3]  =  1.0
	[4, 3]  =  7.0
	[1, 4]  =  6.0
	[2, 4]  =  2.0
	[3, 4]  =  7.0
	
full(gr)

4x4 Array{Float64,2}:
 0.0  3.0  3.0  6.0
 3.0  0.0  1.0  2.0
 3.0  1.0  0.0  7.0
 6.0  2.0  7.0  0.0
~~~

To see the underlying data structure, use `fieldnames`.

~~~julia
fieldnames(gr)

5-element Array{Symbol,1}:
 :m     
 :n     
 :colptr
 :rowval
 :nzval 
~~~

`m` and `n` are the dimensions of the matrix.
The entries of the matrix are stored in nzval.
colptr[i] is the index in nzval of the first nonzero entry
in column i.  rowval tells you which rows in each column are nonzero.
The indices of the nonzero entries in column i are stored in 
rowval[colptr[i]] through rowval[colptr[i+1]-1].

~~~julia
gr.colptr 

5-element Array{Int64,1}:
  1
  4
  7
 10
 13
 
 [gr.rowval gr.nzval]
 
 12x2 Array{Float64,2}:
 2.0  3.0
 3.0  3.0
 4.0  6.0
 1.0  3.0
 3.0  1.0
 4.0  2.0
 1.0  3.0
 2.0  1.0
 4.0  7.0
 1.0  6.0
 2.0  2.0
 3.0  7.0
~~~


