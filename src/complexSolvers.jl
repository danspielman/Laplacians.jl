### samplingSolver
include("revampedLinkedListFloatStorage.jl")
export llsInit
export llsAdd
export llsPurge

include("sqLinOpWrapper.jl")
export SqLinOp

include("samplingSolver.jl")
export samplingParams
export samplingSDDSolver
export samplingLapSolver

### KMP
include("KMPSolver.jl")
export KMPSDDSolver
export KMPLapSolver

### hybridSolver
include("hybridSolver.jl")
export hybridSDDSolver
export hybridLapSolver
export hybridParams

### useful misc
include("condNumber.jl")
export condNumber

include("powerIteration.jl")
export powerIteration

# SDDSolvers take in a SDD matrix - LapSolvers take in an adjacency matrix
SDDSolvers = [augTreeSolver, KMPSDDSolver, hybridSDDSolver, samplingSDDSolver, AMGSolver]
LapSolvers = [augTreeLapSolver, KMPLapSolver, hybridLapSolver, samplingLapSolver, AMGLapSolver]