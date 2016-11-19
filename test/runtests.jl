using Base.Test
using Laplacians

# generate chimeric graphs, and test run many routines on them

function isTree(gr::SparseMatrixCSC)
    isConnected(gr) && (nnz(gr) == 2*(gr.n-1))
end



function testSolvers(a)

    n = a.n
    excess = zeros(n); excess[1] = excess[n] = 0.1;
    la = lap(a);
    sdd = la + spdiagm(excess);
    b = rand(n); b = b - mean(b);

    for solver in SDDSolvers
        f = solver(sdd, tol=1e-6, maxtime=1);
        x = f(b);
        @test norm(sdd*x - b)/norm(b) <= 1e-5
    end

    for solver in LapSolvers
        f = solver(a, tol=1e-6, maxtime=1);
        x = f(b);
        @test norm(la*x - b)/norm(b) <= 1e-5
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
end




include("testPCG.jl")
