# for testing code that is not exported or included in Laplacians.jl

# Test Lex by running code from its notebook

include("$(Pkg.dir("Laplacians"))/src/lex.jl")

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

