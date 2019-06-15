

@deprecate ringGraph(n) ring_graph(n)
@deprecate pathGraph(n) path_graph(n)
@deprecate completeGraph(n) complete_graph(n)
@deprecate generalizedRing(n) generalized_ring(n)
@deprecate randGenRing(n::Int64, k::Integer; verbose=false) rand_gen_ring(n, k; verbose=verbose)
@deprecate hyperCube(n::Int64) hypercube(n)
@deprecate completeBinaryTree(n::Int64) complete_binary_tree(n)

@deprecate productGraph(a0::SparseMatrixCSC, a1::SparseMatrixCSC) product_graph(a0, a1)
#@deprecate wGrid2(n::Integer; weightGen::Function=rand) wgrid2(n, weightGen)
@deprecate randMatching(n::Integer) rand_matching(n)
@deprecate randRegular(n::Integer, k::Integer) rand_regular(n, k)
@deprecate grownGraph(n::Integer, k::Integer) grown_graph(n, k)
@deprecate grownGraphD(n::Integer, k::Integer) grown_graph_d(n, k)
@deprecate prefAttach(n::Integer, k::Integer, p::Float64) pref_attach(n, k, p)

@deprecate randWeight(a) rand_weight(a)
@deprecate randWeightSub(a) rand_weight_sub(a)
@deprecate wtedChimera(n::Integer, k::Integer; verbose=false) wted_chimera(n, k, verbose=verbose)
@deprecate wtedChimera(n::Integer; verbose=false) wted_chimera(n, verbose=verbose)
@deprecate semiWtedChimera(n::Integer; verbose=false, prefix="") semi_wted_chimera(n, verbose=verbose, prefix=prefix)
@deprecate semiWtedChimera(n::Integer, k::Integer; verbose=false, prefix="") semiwted_chimera(n, k, verbose=verbose, prefix=prefix)

@deprecate thickenOnce(a) thicken_once(a)

@deprecate twoLift(a, flip::AbstractArray{Bool,1}) two_lift(a, flip)
@deprecate twoLift(a) two_lift(a)
@deprecate twoLift(a, k::Integer) two_lift(a, k)

@deprecate plotGraph(gr,x,y;color=[0,0,1],dots=true,setaxis=true,number=false) plot_graph(gr,x,y;color=color,dots=dots,setaxis=setaxis,number=number)
@deprecate plotGraph(gr,x,y,z;color=[0,0,1],dots=true,setaxis=true,number=false) plot_graph(gr,x,y,z;color=color,dots=dots,setaxis=setaxis,number=number)

@deprecate approxCholLap(a; args...) approxchol_lap(a; args...)
@deprecate approxCholSDDM(a; args...) approxchol_sddm(a; args...)
@deprecate cholLap(a; args...) chol_lap(a; args...)
@deprecate cholSDDM(a; args...) chol_sddm(a; args...)

@deprecate spectralDrawing(a) spectral_drawing(a)
@deprecate spectralCoords(a) spectral_coords(a)
@deprecate compStretches(tree, mat) comp_stretches(tree, mat)

@deprecate readIJ(fn) read_graph(fn)
@deprecate readIJ(fn, sep) read_graph(fn)
@deprecate readIJV(fn) read_graph(fn)
