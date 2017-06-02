

# testing the interface routines

n = 100;
a = wtedChimera(100,1);
b = randn(n);
b = b - mean(b);
la = Laplacians.forceLap(a);
sdd = copy(la)
sdd[1,1] += 1;
its = [0]

f = Laplacians.lapWrapConnected(cgSolver,a, tol=1e-2, verbose=true)
@test norm(la*f(b)-b)/norm(b) < 1e-1
@test norm(la*f(b,pcgIts=its,tol=1e-3,verbose=false)-b)/norm(b) < 1e-3

norm(la*f(b,pcgIts=its)-b)/norm(b)
norm(la*f(b,pcgIts=its,verbose=true,tol=1e-6)-b)/norm(b)

conSolve = Laplacians.lapWrapConnected(cgSolver)

fa = conSolve(a)

@test norm(la*fa(b,verbose=true) - b) < 1e-2

a2 = disjoin(chimera(100,2),wtedChimera(200,3))
la2 = lap(a2)
n = size(a2,1)
b2 = randn(n);
b2[1:100] = b2[1:100] - mean(b2[1:100]);
b2[101:300] = b2[101:300] - mean(b2[101:300]);

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
    sddm = la + spdiagm(excess);
    b = rand(n); b = b - mean(b);

    its = Int[0]



    for solver in SDDMSolvers
        f = solver(sddm, tol=1e-6, maxtime=maxtime,verbose=true);
        x = f(b,tol=1e-6,maxits=maxits,verbose=true,pcgIts=its);
        @test ((its[1] == maxits) | (norm(sddm*x - b)/norm(b) <= 1e-1))
    end

    for solver in LapSolvers
        f = solver(a, tol=1e-6, maxtime=maxtime,verbose=true);
        x = f(b,tol=1e-6,maxits=maxits,verbose=true,pcgIts=its);
        @test ((its[1] == maxits) | (norm(la*x - b)/norm(b) <= 1e-1))
    end

end



n = 1000
tol = 1e-11

for i in 1:50
    println("wtedChimera($n, $i)")
    if isodd(i)
        gr = wtedChimera(n,i,verbose=true)
    else
        gr = wtedChimera(n,i)
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
    t = akpwU(gru)
end

a = wtedChimera(1234,1);
for i in 2:5
    a = disjoin(a,wtedChimera(1234,i))
end
n = size(a,1)
b = randn(n)
b = b - mean(b)
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
for solver in [augTreeLap, KMPLapSolver, samplingLapSolver, cgLapSolver, ee1, ee2, ee3]
    f = solver(a, tol=1e-6, maxtime=5);
    x = f(b);
    @test norm(la*x - b)/norm(b) <= 1e-1
end



n = 10
a = rand(n,n)
a = a + a'
a = a - diagm(diag(a))
a = sparse(a)
println("rand complete graph")
testSolvers(a,maxtime=10)



n = 20000
i = 1
println("wtedChimera($n, $i)")
gr = wtedChimera(n,i)
testSolvers(gr,maxtime=10)

n = 100000
i = 1
b = randn(n)
b = b - mean(b)
println("wtedChimera($n, $i)")
gr = wtedChimera(n,i)
f = KMPLapSolver(gr,tol=1e-1,verbose=true)
x = f(b)



# Need to add tests for disconnected graphs.
# Not all solvers work with them yet

n = 2000
gr1 = wtedChimera(n,1)
gr2 = wtedChimera(n,2)
gr = disjoin(gr1,gr2)

b1 = rand(n); b1 = b1 - mean(b1)
b2 = rand(n); b2 = b2 - mean(b2)
b = [b1;b2]
f = KMPLapSolver(gr,tol=1e-6)
x = f(b)
la = lap(gr)
@test norm(la*x-b)/norm(b) < 1e-3


# Testing code inside Sampling Solver

a = wtedChimera(100,1)
t = akpw(a)
la = lap(a)
lt = lap(t)
la[1,1] += 1
lt[1,1] += 1
condNumber(la,lt)

sp = Laplacians.defaultSamplingParams

sp.verboseSS = true


n = 100
a = wtedChimera(n,1)
b = randn(n)
b = b - mean(b)
la = lap(a)

f = samplingLapSolver(a,tol=1e-6)
x = f(b)

params = deepcopy(Laplacians.defaultSamplingParams)
params.verboseSS = true
params.returnCN = true
params.fixTree = true
f = samplingLapSolver(a,tol=1e-6,params=params)
x = f(b)





# testing compare_solvers

a = chimera(200,4);

la = lap(a)
sdd = copy(la)
sdd[1,1] += 1
b = randn(a.n)
b = b - mean(b)
@time fn = augTreeLap(a,params=AugTreeParams())
@time fold = augTreeLap(a,params=AugTreeParamsOld())
@time gn = augTreeSddm(sdd)
@time gold = augTreeSddm(sdd,params=AugTreeParamsOld())

lnew(a; kwargs...) = augTreeLap(a; params=AugTreeParams(), kwargs...)
lold(a; kwargs...) = augTreeLap(a; params=AugTreeParamsOld(), kwargs...)
solvers = [SolverTest(lnew,"new") SolverTest(lold,"old")]

dic = Dict()
x = speedTestLapSolvers(solvers, dic, a, b, tol=1e-2)
