type AugTreeHybridParams
    n0::Int64      # the number of edges at which to go direct
    frac::Float64  # add frac*n edges to the tree
    subIts::Float64  # number of its of sub solver
    subTol::Float64  # tol to spec for sub solver
    verbose::Bool
end


