
function isTree(gr::SparseMatrixCSC)
    isConnected(gr) && (nnz(gr) == 2*(gr.n-1))
end

# generate chimeric graphs, and test run many routines on them

function testSolvers(a;maxtime=1)

    n = a.n
    excess = zeros(n); excess[1] = excess[n] = 0.1;
    la = lap(a);
    sdd = la + spdiagm(excess);
    b = rand(n); b = b - mean(b);

    for solver in SDDSolvers
        f = solver(sdd, tol=1e-6, maxtime=maxtime);
        x = f(b);
        @test norm(sdd*x - b)/norm(b) <= 1e-1
    end

    for solver in LapSolvers
        f = solver(a, tol=1e-6, maxtime=maxtime);
        x = f(b);
        @test norm(la*x - b)/norm(b) <= 1e-1
    end

end



n = 1000
tol = 1e-11

for i in 1:100
    println("wtedChimera($n, $i)")
    gr = wtedChimera(n,i)
    @test gr.n == n
    @test isConnected(gr)
    @test isTree(kruskal(gr))
    @test isTree(randishKruskal(gr))
    @test isTree(randishPrim(gr))
    @test abs(sum(kruskal(gr)) - sum(prim(gr))) < tol
    @test sum(kruskal(gr,kind=:min)) <= sum(kruskal(gr,kind=:max))

    testSolvers(gr)

    gru = unweight(gr)
    t = akpwU(gru)
end

n = 20000
i = 1
println("wtedChimera($n, $i)")
gr = wtedChimera(n,i)
testSolvers(gr,maxtime=10)

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

