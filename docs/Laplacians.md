[TOC]

# yinsGraph

yinsGraph is a package that I (Dan) am writing to explore and manipulate graphs in Julia.  The graphs are represented as sparse matrices.  The particular class in Julia is called a SparseMatrixCSC.  The reasons for this are:

  * They are fast, and
  * We want to do linear algebra with them, so matrices help.
  
You can probably learn more about the CSC (Compressed Sparse Column) format by googling it.  

So far, speed tests of code that I've written for connected components, shorest paths, and minimum spanning trees have been as fast or faster than the previous routines I could call from Matlab.  

  

## To install yinsGraph
you will need a number of packages.
You install these like
```
Pkg.add("PyCall")
Pkg.add("PyPlot")
Pkg.add("DataStructures")
```
I also recommend the Optim package.

I think you need to install matplotlib in python before PyPlot.  
Look at this page for more information: https://github.com/stevengj/PyPlot.jl

I'm not sure if there are any others.  If you find that there are, please list them above.

## To use yinsGraph

Examples of how to do many things in yinsGraph may be found in the IJulia notebooks.  These have the extensions .ipynb.  When they look nice, I think it makes sense to convert them to .html.  

Right now, the notebooks worth looking at are:

* [yinsGraph](yinsGraph.html) - usage, demo, and speed tests
* [Solvers](Solvers.html) - code for solving equations.  How to use direct methods, conjugate gradient, and a preconditioned augmented spanning tree solver.

* [ ] The implementation of CG in IterativeSolvers sort of sucks, as I now see in the tests.  It is allocating way to much memory.  It should be fixed either by devectorizing (see Julia Performance Tips), or by using BLAS routines.

(I suggest that you open the html in your browser)

### Graph generators:
~~~julia
 readIJ(filename::String)
 readIJV(filename::String)
 writeIJV(filename::String, mat)
 ringGraph(n::Int64)
 generalizedRing(n::Int64, gens)
 randMatching(n::Int64)
 randRegular(n::Int64, k::Int64)
 grownGraph(n::Int64, k::Int64)
 grownGraphD(n::Int64, k::Int64)
 prefAttach(n::Int64, k::Int64, p::Float64)
 hyperCube(d::Int64)
 completeBinaryTree(n::Int64)
 grid2(n::Int64)
 grid2(n::Int64, m::Int64; isotropy=1)
 grid2coords(n::Int64, m::Int64)
~~~

* [ ] The types in the arguments of those should probably be more flexible/general.

For example, to generate a 4-by-5 grid, you type

~~~julia
graph = grid2(4,5)
~~~

### Operations on Graphs:

* `shortIntGraph`  for converting the index type of a graph to an Int32.  
* `lap`  to produce the laplacian of a graph
  * [ ] Maybe this should grab the upper triangular part, and symmetrize first. 	
* `unweight` - change all the weights to 1
* `mapweight{Tval,Tind}(a::SparseMatrixCSC{Tval,Tind},f)`  to apply the function f to the weight of every edge.
* `uniformWeight`  an example of mapweight.  It ignores the weight, and maps every weight to a random in [0,1]
* `productGraph(a0::SparseMatrixCSC, a1::SparseMatrixCSC)` the cartesian product.  Given two paths it makes a grid.
* `edgeVertexMat(mat::SparseMatrixCSC)`  signed edge vertex matrix
* `subsampleEdges(a::SparseMatrixCSC{Float64,Int64}, p::Float64)`
  produce a new graph that keeps each edge with probability p.
* `twoLift(a, k)` create a 2-lift of a with k flipped edges.  If k is unspecified, this generates a random 2-lift. 
* `joinGraphs(a, b, k)` create a disjoint union of a and b, and add k random edges between them
* `plotGraph(gr,x,y,color=[0,0,1];dots=true,setaxis=true,number=false)`
* `spectralDrawing(graph)`

### Fundamental Graph Algorithms:

*  `components` computes connected components, returns as a vector
*  `vecToComps` turns into an array with a list of vertices in each component
*  `shortestPaths(mat, start)`  returns an array of distances,
    and pointers to the node closest (parent array)
*  `kruskal(mat; kind=:min)`  to get a max tree, use `kind = :max`
    returns it as a sparse matrix.
 
### Solving Linear equations:

We have implemented Conjugate Gradient (cg) and the Preconditioned Conjugate Gradient (pcg).  These implementations use BLAS when they can, and a slower routine for data types like BigFloat.

To learn more, read [solvers.md](solvers.md).


## To develop yinsGraph

Just go for it.
Don't worry about writing fast code at first.
Just get it to work.
We can speed it up later.
The yinsGraph.ipynb notebook contains some examples of speed tests.
Within some of the files, I am keeping old, unoptimized versions of code around for comparison (and for satisfaction).  I will give them the name "XSlow"

I think that each file should contain a manifest up top listing the functions and types that it provides.  They should be divided up into those that are for internal use only, and those that should be exported.  Old code that didn't work well, but which you want to keep for reference should go at the end.

### Using sparse matrices as graphs
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

### Parametric Types

A sparse matrix has two types associated with it: the types of its indices (some sort of integer) and the types of its values (some sort of number).  Most of the code has been written so that once these types are fixed, the type of everything else in the function has been too.  This is accomplished by putting curly braces after a function name, with the names of the types that we want to use in the braces.  For example,

~~~julia
shortestPaths{Tv,Ti}(mat::SparseMatrixCSC{Tv,Ti}, start::Ti)
~~~

`Tv`, sometimes written `Tval` denotes the types of the values, and `Ti` or `Tind` denotes the types of the indices.  This function will only be called if the node from which we compute the shortest paths, `start` is of type `Ti`.  Inside the code, whenever we write something like `pArray = zeros(Ti,n)`, it creates an array of zeros of type Ti.  Using these parameteric types is *much* faster than leaving the types unfixed.

### Data structures:

* `IntHeap` a heap that stores small integers (like indices of nodes in a graph) and that makes deletion fast.  Was much faster than using Julia's more general heap.

### Interface issue:
There are many different sorts of things that our code could be passing around.  For example, kruskal returns a graph as a sparse matrix.  But, we could use a format that is more specialized for trees, like the RootedTree type.  At some point, when we optimize code, we will need to figure out the right interfaces between routines.  For example, some routines symmetrize at the end.  This is slow, and should be skipped if not necessary.  It also doubles storage.

### Writing tests:
I haven't written any yet.  I'll admit that I'm using the notebooks as tests.  If I can run all the cells, then it's all good.


