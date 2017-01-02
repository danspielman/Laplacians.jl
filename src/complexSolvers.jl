### samplingSolver
include("revampedLinkedListFloatStorage.jl")
export llsInit
export llsAdd
export llsPurge

include("sqLinOpWrapper.jl")
export SqLinOp

include("samplingSolver.jl")
export samplingSDDMSolver
export samplingLapSolver
export samplingParams

### KMP
# included in Laplacians.jl

### hybridSolver
include("hybridSolver.jl")
export hybridSDDMSolver
export hybridLapSolver
export hybridParams

### useful misc
include("condNumber.jl")
export condNumber

include("powerIteration.jl")
export powerIteration

"""A list containing SDDM linear system solvers. They take in a SDDM matrix plus tol, maxits and maxtime parameters."""
SDDMSolvers = [augTreeSolver, KMPSDDMSolver, samplingSDDMSolver, AMGSolver]

"""A list containing Laplacian linear system solvers. They take in an adjacency matrix plus tol, maxits and maxtime parameters."""
LapSolvers = [augTreeLapSolver, KMPLapSolver, samplingLapSolver, AMGLapSolver, cgLapSolver]
