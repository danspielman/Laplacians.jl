

# testing the interface routines

n = 100;
a = wted_chimera(100,1);
b = randn(n);
b = b .- mean(b);
la = Laplacians.forceLap(a);
sdd = copy(la)
sdd[1,1] += 1;
its = [0]

f = Laplacians.lapWrapConnected(cgSolver,a, tol=1e-2, verbose=true)
@test norm(la*f(b)-b)/norm(b) < 1e-1
@test norm(la*f(b,pcgIts=its,tol=1e-3,verbose=false)-b)/norm(b) < 1e-2

f(zeros(n))

norm(la*f(b,pcgIts=its)-b)/norm(b)
norm(la*f(b,pcgIts=its,verbose=true,tol=1e-6)-b)/norm(b)

conSolve = Laplacians.lapWrapConnected(cgSolver)

fa = conSolve(a)

@test norm(la*fa(b,verbose=true) - b) < 1e-2

a2 = disjoin(chimera(100,2),wted_chimera(200,3))
la2 = lap(a2)
n = size(a2,1)
b2 = randn(n);
b2[1:100] = b2[1:100] .- mean(b2[1:100]);
b2[101:300] = b2[101:300] .- mean(b2[101:300]);

f = Laplacians.lapWrapComponents(conSolve, a2)

x = f(b2,verbose=true)
norm(la2*x-b2)/norm(b2)

f = Laplacians.lapWrapComponents(Laplacians.cgLapSolver, a2)

x1 = f(b2,verbose=true)
norm(la2*x1-b2)/norm(b2)

sum(x1[101:300])

f = Laplacians.cgLapSolver(a)
x = f(b,verbose=true)
norm(la*x-b)/norm(b)

f0 = Laplacians.lapWrapComponents(conSolve, a)
@test norm(la*f0(b,tol=1e-4,maxits=200,verbose=true)-b)/norm(b) < 1e-2

f0 = Laplacians.lapWrapComponents(Laplacians.cgLapSolver, a)
@test norm(la*f0(b,tol=1e-4,maxits=200,verbose=true)-b)/norm(b) < 1e-2


f = Laplacians.lapWrapComponents(Laplacians.lapWrapConnected(cholSDDM),a)
@test norm(la*f(b)-b)/norm(b) < 1e-8

solver = Laplacians.lapWrapComponents(Laplacians.lapWrapConnected(Laplacians.cholSDDM))
f = solver(a,verbose = true)
norm(la*f(b,tol=1e-3)-b)/norm(b)

its2 = [0]
f = Laplacians.lapWrapSDDM(Laplacians.cgSolver,a,tol=1e-2,pcgIts=its)
@test norm(la*f(b,verbose=true,pcgIts = its2, tol=1e-3)-b) / norm(b) < 1e-2
f = Laplacians.lapWrapSDDM(Laplacians.cholSDDM,a)
@test norm(la*f(b)-b) / norm(b) < 1e-2
f = Laplacians.cholLap(a)
@test norm(la*f(b)-b) / norm(b) < 1e-2
fs = Laplacians.lapWrapSDDM(cgSolver)
f = fs(a, tol=1e-2, verbose=true)
x = f(b, tol=1e-6);




mats = []
rhss = []
solver = Laplacians.wrapCapture(approxCholLap, mats, rhss)
a = chimera(10)
as = shortIntGraph(a)
f = solver(a);
fs = solver(as);
size(mats[1])
b = randn(10)
x = f(b);
xs = fs(b);

solver = Laplacians.approxCholLapChol(a,verbose=true)
x = solver(b);

# testing approxChol internals
a = randRegular(20,3)
llp = Laplacians.LLmatp(a)
Laplacians.print_ll_col(llp,1)
llo = Laplacians.LLMatOrd(a);
Laplacians.print_ll_col(llo,1)

a = wted_chimera(1111,1);
llp = Laplacians.LLmatp(a)
Random.seed!(1)
ldli = Laplacians.approxChol(llp);
@test Laplacians.condNumber(a, ldli, verbose=true) < 20


# testing by repitition

function isTree(gr::SparseMatrixCSC)
    isConnected(gr) && (nnz(gr) == 2*(gr.n-1))
end

# generate chimeric graphs, and test run many routines on them

function testSolvers(a;maxtime=5)

    maxits = 200

    n = a.n
    excess = zeros(n); excess[1] = excess[n] = 0.1;
    la = lap(a);
    sddm = la + sparse(Diagonal(excess));
    b = rand(n); b = b .- mean(b);

    its = Int[0]



    for solver in SDDMSolvers
        println(solver)
        f = solver(sddm, tol=1e-6, maxtime=maxtime,verbose=true);
        x = f(b,tol=1e-6,maxits=maxits,verbose=true,pcgIts=its);
        @test ((its[1] == maxits) | (norm(sddm*x - b)/norm(b) <= 1e-1))
    end

    for solver in LapSolvers
        println(solver)
        f = solver(a, tol=1e-6, maxtime=maxtime,verbose=true);
        x = f(b,tol=1e-6,maxits=maxits,verbose=true,pcgIts=its);
        @test ((its[1] == maxits) | (norm(la*x - b)/norm(b) <= 1e-1))
    end

end



n = 1000
tol = 1e-11

for i in 1:50
    println()
    println("wted_chimera($n, $i)")
    if isodd(i)
        gr = wted_chimera(n,i,verbose=true, ver=L.V06)
    else
        gr = wted_chimera(n,i, ver=L.V06)
    end
    @test gr.n == n
    @test isConnected(gr)
    @test isTree(kruskal(gr))
    @test isTree(randishKruskal(gr))
    @test isTree(randishPrim(gr))
    @test abs(sum(kruskal(gr)) - sum(prim(gr))) < tol
    @test sum(kruskal(gr,kind=:min)) <= sum(kruskal(gr,kind=:max))

    nnzL, flops = ask_cholmod(lap(gr))
    pe = cholmod_perm(lap(gr))

    testSolvers(gr)

    gru = unweight(gr)
    _ = akpwU(gru)

end

for i in 1:50
    wted_chimera(n, i, verbose=true)
end

a = wted_chimera(1234,1);
for i in 2:5
    global a = disjoin(a,wted_chimera(1234,i))
end
n = size(a,1)
b = randn(n)
b = b .- mean(b)
la = lap(a)

ee1 = function(a; verbose=false, args...)
    approxCholLap(a; params=ApproxCholParams(:deg), verbose=verbose, args...)
end
ee2 = function(a; verbose=false, args...)
    approxCholLap(a; params=ApproxCholParams(:wdeg), verbose=verbose, args...)
end
ee3 = function(a; verbose=false, args...)
    approxCholLap(a; params=ApproxCholParams(:given), verbose=verbose, args...)
end
for sol in [augTreeLap, KMPLapSolver, cgLapSolver, ee1, ee2, ee3]
    fx1 = sol(a, tol=1e-6, maxtime=5);
    xx1 = fx1(b);
    @test norm(la*xx1 - b)/norm(b) <= 1e-1
end



n = 10
a = rand(Float64,n,n)
a = a + a'
a = a - Diagonal(diag(a))
a = sparse(a)
println("rand complete graph")
testSolvers(a,maxtime=10)

n = 100
a = rand(Float64, n,n)
a = a + a'
a = a - Diagonal(diag(a))
a = sparse(a)
println("rand complete graph")
testSolvers(a,maxtime=10)


n = 20000
i = 1
println("wted_chimera($n, $i)")
gr = wted_chimera(n,i)
testSolvers(gr,maxtime=10)

n = 100000
i = 1
b = randn(n)
b = b .- mean(b)
println("wted_chimera($n, $i)")
gr = wted_chimera(n,i)
f = KMPLapSolver(gr,tol=1e-1,verbose=true)
x = f(b)



# Need to add tests for disconnected graphs.
# Not all solvers work with them yet

n = 2000
gr1 = wted_chimera(n,1)
gr2 = wted_chimera(n,2)
gr = disjoin(gr1,gr2)

b1 = rand(n); b1 = b1 .- mean(b1)
b2 = rand(n); b2 = b2 .- mean(b2)
b = [b1;b2]
f = KMPLapSolver(gr,tol=1e-6)
x = f(b)
la = lap(gr)
@test norm(la*x-b)/norm(b) < 1e-3


# Testing code inside Sampling Solver

a = wted_chimera(100,1)
t = akpw(a)
la = lap(a)
lt = lap(t)
la[1,1] += 1
lt[1,1] += 1
condNumber(la,lt)
condNumber(la,lt,tol=1)

gOp = Laplacians.SqLinOp(true,1.0,100,x->la*x)
R = powerIteration(gOp, tol=1e-8, maxit=2, verbose=true)
R = powerIteration(gOp, tol=1e-1, maxit=100, verbose=true)


#=
sp = Laplacians.defaultSamplingParams

sp.verboseSS = true


n = 100
a = wted_chimera(n,1)
b = randn(n)
b = b .- mean(b)
la = lap(a)

f = samplingLapSolver(a,tol=1e-6)
x = f(b)

params = deepcopy(Laplacians.defaultSamplingParams)
params.verboseSS = true
params.returnCN = true
params.fixTree = true
f = samplingLapSolver(a,tol=1e-6,params=params)
x = f(b)
=#




# testing compare_solvers

a = chimera(200,4,ver=L.V06);

la = lap(a)
sdd = copy(la)
sdd[1,1] += 1
b = randn(a.n)
b = b .- mean(b)
@time fn = augTreeLap(a,params=AugTreeParams())
@time fold = augTreeLap(a,params=AugTreeParamsOld())
@time gn = augTreeSddm(sdd)
@time gold = augTreeSddm(sdd,params=AugTreeParamsOld())

lnew(a; kwargs...) = augTreeLap(a; params=AugTreeParams(), kwargs...)
lold(a; kwargs...) = augTreeLap(a; params=AugTreeParamsOld(), kwargs...)
solvers = [SolverTest(lnew,"new") SolverTest(lold,"old")]

dic = Dict()
x = speedTestLapSolvers(solvers, dic, a, b, tol=1e-2)
