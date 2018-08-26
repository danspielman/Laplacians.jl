### samplingSolver
include("revampedLinkedListFloatStorage.jl")
export llsInit
export llsAdd
export llsPurge

include("sqLinOpWrapper.jl")
export SqLinOp

#=
include("samplingSolver.jl")
export samplingSDDMSolver
export samplingLapSolver
export samplingParams
=#

### KMP
# included in Laplacians.jl

### hybridSolver -- has been removed from master, and now lives in own branch

### useful misc
include("condNumber.jl")
export condNumber

include("powerIteration.jl")
export powerIteration

"""A list containing SDDM linear system solvers. They take in a SDDM matrix plus tol, maxits and maxtime parameters."""
SDDMSolvers = [augTreeSddm, KMPSDDMSolver, approxchol_sddm]

"""A list containing Laplacian linear system solvers. They take in an adjacency matrix plus tol, maxits and maxtime parameters."""
LapSolvers = [approxchol_lap, augTreeLap, KMPLapSolver, cgLapSolver]
