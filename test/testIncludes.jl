# for testing code that is not exported or included in Laplacians.jl



#============================================================#
# Test compare_solvers code

#=
include("$(Pkg.dir("Laplacians"))/src/compare_solvers.jl")

lnew(a; kwargs...) = augTreeLap(a; params=AugTreeParams(), kwargs...)
lold(a; kwargs...) = augTreeLap(a; params=AugTreeParamsOld(), kwargs...)
solvers = [SolverTest(lnew,"new") SolverTest(lold,"old")]

dic = Dict()
x = speedTestLapSolvers(solvers, dic, a, b, tol=1e-2)
=#

#============================================================#
# Test Lex by running code from its notebook

include("$(Pkg.dir("Laplacians"))/src/lex.jl")


setLexDebugFlag(true)

n = 10
Pn = pathGraph(n)

isTerm = zeros(Bool, n)
isTerm[1] = true
isTerm[n] = true

initVal = zeros(n)
initVal[n] = 1.0

# inf-minimizer
infMinVolt = CompInfMin(Pn, isTerm, initVal)
println(infMinVolt)
println(MaxEdgeGrad(Pn, infMinVolt))

# lex-minimizer
lexMinVolt = CompLexMin(Pn, isTerm, initVal)

println(lexMinVolt)

n = 100
G = chimera(n,1)
isTerm = zeros(Bool, n)
# arbitrary terminal values
isTerm[1] = true
isTerm[5] = true
isTerm[11] = true
isTerm[18] = true

initVal = zeros(Float64, n)
initVal[1] = 0.0
initVal[5] = 13
initVal[11] = 7
initVal[18] = 11

infMinVolt = CompInfMin(G, isTerm, initVal)
println(infMinVolt)
println(MaxEdgeGrad(G, infMinVolt))

lexMinVolt = simIterLex(500, G, isTerm, initVal)
println(lexMinVolt)
println(MaxEdgeGrad(G, lexMinVolt))
println(checkLex(G, isTerm, initVal, lexMinVolt))

simIterLexTest()

#============================================================#
# test isotonicIPM

include("$(Pkg.dir("Laplacians"))/src/isotonicIPM.jl")

n = 100;
A = grid2(10,10)
A = triu(A);
v = randn(n)+collect(1:n)/2;
isoTestValues = collect(1:n);
newOrder = randperm(n)
permMat = (speye(n))[1:n,newOrder];
A = permMat*A*permMat';
v = permMat*v;
@time (x,acc,itercount) = isotonicIPM(A,v);

