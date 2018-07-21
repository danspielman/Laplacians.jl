

@deprecate ringGraph(n) ring_graph(n)
@deprecate pathGraph(n) path_graph(n)
@deprecate completeGraph(n) complete_graph(n)
@deprecate generalizedRing(n) generalized_ring(n)
@deprecate randGenRing(n::Int64, k::Integer; verbose=false) rand_gen_ring(n, k; verbose=verbose) 
@deprecate hyperCube(n::Int64) hypercube(n)
@deprecate completeBinaryTree(n::Int64) complete_binary_tree(n)

@deprecate productGraph(a0::SparseMatrixCSC, a1::SparseMatrixCSC) product_graph(a0, a1)
@deprecate wGrid2(n::Integer; weightGen::Function=rand) wgrid2(n, weightGen)
@deprecate randMatching(n::Integer) rand_matching(n)
@deprecate randRegular(n::Integer, k::Integer) rand_regular(n, k)
@deprecate grownGraph(n::Integer, k::Integer) grown_graph(n, k)
@deprecate grownGraphD(n::Integer, k::Integer) grown_graph_d(n, k)
@deprecate prefAttach(n::Integer, k::Integer, p::Float64) pref_attach(n, k, p)
