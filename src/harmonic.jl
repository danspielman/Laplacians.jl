"""
    x = harmonic_interp(adj_mat, S, vals; tol=1e-6)

Interpolates a function on a graph, given by its adjacency matrix,
by minizing the Laplacian quadratic form subject to the boundary
conditions that `x[S[i]] = vals[i]` for `i` in `S`.

This is the algorithm sometimes known as Label Propagation,
or Semi-Supervised Learning on Graphs.  The idea comes from the paper
"Semi-Supervised Learning Using Gaussian Fields and Harmonic Functions"
by Zhu, Gharamani, and Lafferty from ICML 2003.

This version might fail for disconnected graphs.
You can check if a graph is connected with `isConnected(adj_mat)`.
"""
function harmonic_interp(adj_mat, S::Vector, vals::Vector; tol=1e-6)
    n = size(adj_mat,1)
    b = zeros(n)
    b[S] = vals

    inds = ones(Bool,n)
    inds[S] .= false
    la = lap(adj_mat)
    la_sub = la[inds,inds]
    b_sub = (-la*b)[inds]

    f = approxchol_sddm(la_sub; tol=tol)
    x_sub = f(b_sub)

    x = copy(b)
    x[inds] = x_sub
    return x
end
