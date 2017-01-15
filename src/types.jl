

"""
Parameters for the sampling solver.
"""
type SimpleSamplerParams{Tv,Ti}
    startingSize::Ti    # the initial size of the linked list storage structure
    blockSize::Ti       # the size of each consecutive block of memory assigned to a certain element

    # debug parameters, we can get rid of them later
    verboseSS::Bool
    verbose::Bool
    returnCN::Bool
    CNTol::Tv
    order::Symbol    # options :min, :approx, :tree
end

type AugTreeHybridParams
    n0::Int64      # the number of edges at which to go direct
    frac::Float64  # add frac*n edges to the tree
    subIts::Float64  # number of its of sub solver
    subTol::Float64  # tol to spec for sub solver
    verbose::Bool
    sampParams::SimpleSamplerParams
end
