

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
@assert norm(la*f(b)-b)/norm(b) < 1e-1
@assert norm(la*f(b,pcgIts=its,tol=1e-3,verbose=false)-b)/norm(b) < 1e-3

norm(la*f(b,pcgIts=its)-b)/norm(b)
norm(la*f(b,pcgIts=its,verbose=true,tol=1e-6)-b)/norm(b)

conSolve = Laplacians.lapWrapConnected(cgSolver)

fa = conSolve(a)

@assert norm(la*fa(b,verbose=true) - b) < 1e-2

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
@assert norm(la*f0(b,tol=1e-4,maxits=200,verbose=true)-b)/norm(b) < 1e-2

f0 = Laplacians.lapWrapComponents(Laplacians.cgLapSolver, a)
@assert norm(la*f0(b,tol=1e-4,maxits=200,verbose=true)-b)/norm(b) < 1e-2


f = Laplacians.lapWrapComponents(Laplacians.lapWrapConnected(cholSDDM),a)
@assert norm(la*f(b)-b)/norm(b) < 1e-8

solver = Laplacians.lapWrapComponents(Laplacians.lapWrapConnected(Laplacians.cholSDDM))
f = solver(a,verbose = true)
norm(la*f(b,tol=1e-3)-b)/norm(b)

its2 = [0]
f = Laplacians.lapWrapSDDM(Laplacians.cgSolver,a,tol=1e-2,pcgIts=its)
@assert norm(la*f(b,verbose=true,pcgIts = its2, tol=1e-3)-b) / norm(b) < 1e-2
f = Laplacians.lapWrapSDDM(Laplacians.cholSDDM,a)
@assert norm(la*f(b)-b) / norm(b) < 1e-2
f = Laplacians.cholLap(a)
@assert norm(la*f(b)-b) / norm(b) < 1e-2
fs = Laplacians.lapWrapSDDM(cgSolver)
f = fs(a, tol=1e-2, verbose=true)
x = f(b, tol=1e-6);
             


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

for i in 1:100
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
for solver in [augTreeLap, KMPLapSolver, samplingLapSolver, cgLapSolver]
    f = solver(a, tol=1e-6, maxtime=5);
    x = f(b);
    @test norm(la*x - b)/norm(b) <= 1e-1
end




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

# One test for Hybrid Solver, just to get coverage
n = 100
a = wtedChimera(n,1)
b = randn(n)
b = b - mean(b)
la = lap(a)
f = hybridLapSolver(a,tol=1e-6)
x = f(b)
sddm = la + diagm(rand(n)/1000)
f = hybridSDDMSolver(sddm,tol=1e-6)
x = f(b)

hp = Laplacians.defaultHybridParams
hp.n0=10
f = hybridSDDMSolver(sddm,tol=1e-6,verbose=true)
x = f(b)

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

