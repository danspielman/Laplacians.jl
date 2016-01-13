# treeAlgs
### compStretches
Compute the stretched of every edge in `mat` with respect to the tree `tree`. Returns the answer as a sparse matrix with the same nonzero structure as `mat`. Assumes that `mat` is symmetric. `tree` should be the adjacency matrix of a spanning tree.


```julia
compStretches{Tv,Ti}(t::Laplacians.RootedTree{Tv,Ti}, mat::SparseMatrixCSC{Tv,Ti})
compStretches{Tv,Ti}(tree::SparseMatrixCSC{Tv,Ti}, mat::SparseMatrixCSC{Tv,Ti})
```

treeAlgs.jl:393



