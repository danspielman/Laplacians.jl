# Test every function exported in Laplacians.jl by running it at least once
# Create this by copying over the export commands from Laplacians.jl

n = 101
a = wted_chimera(n,1)

# export symPermuteCSC

rp = randperm(n)
ap = symPermuteCSC(a,rp);
rpi = zeros(Int,n);
rpi[rp] = 1:n

@test sum(abs.(a-symPermuteCSC(ap,rpi))) == 0

# export symTransposeCSC

a2 = triu(a) + (tril(a) .> 0)
@test sum(abs.(symTransposeCSC(a2) - a2')) == 0

# export submatrixCSC

s = collect(1:10)
submatrixCSC(a,s)

  # export deg
  # export nbri
  # export weighti
  # export nbrs
  # export wdeg

x = 0
y = 0
z = 0
w = 0
for i in 1:n
    for j in 1:deg(a,i)
        global x += weighti(a,i,j)
        global y += a[i,nbri(a,i,j)]
    end
    for j in nbrs(a,i)
        global z += a[i,j]
    end
    global w += wdeg(a,i)
end
@test isapprox(x,y)
@test isapprox(x,z)
@test isapprox(x,w)
@test isapprox(x,sum(a))

# export setValue

setValue(a,1,1,0.0)

# export backIndices

b = backIndices(a)

# export flipIndex

b = flipIndex(a)

# export findEntries

u,v,w = findEntries(a)

# export compConductance

compConductance(a,collect(1:10))

# export getVolume

getVolume(a,collect(1:10))

# export getObound

getObound(a,collect(1:10))

  # export readIJ
  # export ringGraph
  # export generalizedRing
  # export randMatching
  # export rand_regular


# export grown_graph

a2 = grown_graph(100,3)

# export grown_graph_d

a2 = grown_graph_d(100,3)

# export pref_attach

a2 = pref_attach(100,3,0.5)
a2 = pref_attach(5,4,0.5)


# export hypercube

a2 = hypercube(3)

# export complete_binary_tree

a2 = complete_binary_tree(7)

# export complete_graph

a2 = complete_graph(7)

# export complete_bipartite_graph

a2 = complete_bipartite_graph(7)

# export path_graph

a2 = path_graph(7)

# export star_graph

a2 = star_graph(7)

# export grid2

a2 = grid2(3)

# export grid2coords

a2 = grid2coords(3)

# export grid3

a3 = grid3(3)

  # export randGenRing
  # export randperm

# export ErdosRenyi

a2 = ErdosRenyi(100,300)

# export ErdosRenyiCluster

a2 = ErdosRenyiCluster(100,4)

# export ErdosRenyiClusterFix

a2 = ErdosRenyiClusterFix(100,4)

# export pureRandomGraph

# export chimera
# export randWeight
# export wted_chimera, semiwted_chimera

for i in 1:5
    semiwted_chimera(10000,i)
end


# export readIJ, read_graph, writeIJV
println("Testing IO")
n = 101
a = wted_chimera(n,1)
writeIJV("tmp.txt",a)
a2 = read_graph("tmp.txt")
@test sum(abs.(a-a2)) == 0

a2 = read_graph("tmp.txt")
@test sum(abs.(a-a2)) == 0

fh = open("tmp.txt","w")
write(fh,"1, 3, 4 \n 2, 3, 2.5 \n")
close(fh)

a1 = read_graph("tmp.txt")

fh = open("tmp.txt","w")
write(fh,"1 3 4 \n 2 3 2.5 \n")
close(fh)

a2 = read_graph("tmp.txt")

@test a1 == a2

rm("tmp.txt")


# export unweight, unweight!

a2 = unweight(a2)
unweight!(a2)

  # export mapweight
  # export uniformWeight, uniformWeight!

a2 = uniformWeight(a2)
uniformWeight!(a2)

  # export edgeVertexMat

b = edgeVertexMat(a2)
@test sum(abs.(b'*b - lap(unweight(a2)))) == 0

  # export line_graph

c = line_graph(path_graph(10))
@test nnz(c) == 16
c = line_graph(ring_graph(10))
@test nnz(c) == 20
c = line_graph(complete_graph(10))
@test nnz(c) == 16*10*9/2

a = wted_chimera(102,2)
b = wtedEdgeVertexMat(a)
@test sum(abs.(b'*b - lap(a))) < 1e-8

  # export power, thicken_once, thicken

  a = power(grid2(10),4)
  a = thicken_once(grid2(10))
  a = thicken(grid2(10),4)

  # export productGraph
  # export generalizedNecklace
  # export subsampleEdges

subsampleEdges(a,0.5)

# export two_lift

two_lift(a)
two_lift(a,3)

  # export joinGraphs, disjoin

  # export plotGraph

  # export shortIntGraph, floatGraph

  # export lap
  # export adj
  # export spectral_coords
  # export spectral_drawing

  # export toUnitVector

# export diagmat

diagmat(a)

  # export components
  # export biggestComp
  # export vecToComps
  # export isConnected

# export shortestPaths, shortestPathTree, pathFromParents

a = wted_chimera(102,1)
shortestPaths(a,1)
shortestPathTree(a,1)

Laplacians.intHeapSort(randn(10))

nh = Laplacians.intHeap(10)
for i in 1:10
    Laplacians.intHeapAdd!(nh, i, rand())
end
Laplacians.intHeapSort(nh)

  # export kruskal, prim


# export RootedTree
# export matToTree

# export matToTreeDepth

a = wted_chimera(101,1)
t = akpw(a)
tr = matToTree(t)
tr, d1 = matToTreeDepth(t);
d2 = Laplacians.treeDepthDFS(t)

  # export tarjanStretch
  # export compDepth
  # export comp_stretches
  # export dfsOrder

Laplacians.bfsOrder(t,1)

t0 = complete_binary_tree(6); t0[1,2] = 0; t0[2,1] = 0; dropzeros!(t0);
@test_throws ErrorException Laplacians.bfsOrder(t0, 1)
@test_throws ErrorException Laplacians.matToTree(t0)
@test_throws ErrorException Laplacians.matToTreeDepth(t0)

  # export cg, cgSolver
  # export pcg, pcgSolver, pcgLapSolver

  # export maxflow

for i in 1:10
  a = wted_chimera(100,i)
  a[90:99,100] .= 2;
  a[100,90:99] .= 2;
  a[1,2:10] .= 2;
  a[2:10,1] .= 2;
  f,c = maxflow(a,1,100)

  @test sum(abs.(f+f')) < 1e-8
  y = f*ones(100)
  y[1] = 0
  y[100] = 0
  @test sum(abs.(y)) < 1e-8

  x = zeros(100)
  x[c] .= 1.0
  t = findall(iszero,x)
  @test abs( sum(a[c,t]) - sum(f[1,:])  ) < 1e-8

  @test maximum(f - a) < 1e-8

end

  # export akpw, akpwU

akpw(wted_chimera(10000,1),ver=2)

  # export prn
  # export apr
  # export localImprove
  # export refineCut
  # export dumb

a = chimera(100, 3);
s = prn(a, [1,2,3], 0.2, 5);
conds = compConductance(a, s)
#println(conds, " ", length(s))
minEpsSigma = getVolume(a, s) / getVolume(a, setdiff(collect(1:max(a.n, a.m)), s));
cut, flow = localImprove(a, s, epsSigma = minEpsSigma);
condcut = compConductance(a, cut)
heur = refineCut(a, cut)
dumbRefineCut(a,collect(1:10))

#println(condcut, " ", length(cut))

  # export randishKruskal, randishPrim

  # export FastSampler, sample, sampleMany


r = rand(10)
s = FastSampler(r)
sample(s)
blockSample(r)

# export fiedler

fiedler(chimera(100))

  # export SolverTest, speedTestLapSolvers

solvers = [SolverTest(approxchol_lap,"ac") SolverTest(augTreeLap,"aug")]

dic = Dict()
n = 1000
a = chimera(n)
b = randn(n)
b = b .- mean(b)
x = speedTestLapSolvers(solvers, dic, a, b, tol=1e-2, verbose=true)

f = Laplacians.augTreeFactor(a, akpw(a));

  # include("conditionNumber.jl")
  # export support, approxQual, conditionNumber

  # include("sparsify.jl")
  # export sparsify

a = wted_chimera(1000,1)
dave = nnz(a)/size(a,1)
a = thicken(a,round(Int,200/dave))
as = sparsify(a,ep=1,JLfac=4);
@test approxQual(a,as,verbose=true) < 2
@test conditionNumber(a,as,tol=1e-4) < 10



  # export johnlind

a = chimera(n,1)
johnlind(a)

  # export toposort, dirEdgeVertexMat

a = wted_chimera(301,1)
S = [1;2;3]
vals = [1.0;2.0;3.0]
x = harmonic_interp(a, S, vals, tol=1e-10)
@test x[S] == vals
b = lap(a)*x
b[S] .= 0
@test sum(abs.(b)) < 1e-6




println("End of testByExport")
