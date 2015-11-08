using Base.Test

# generate chimeric graphs, and test run many routines on them

function isTree(gr::SparseMatrixCSC)
    isConnected(gr) && (nnz(gr) == 2*(gr.n-1))
end

n = 1000

for i in 1:100
    println("wtedChimera($n, $i)")
    gr = wtedChimera(n,i)
    @test gr.n == n
    @test isConnected(gr)
    @test isTree(kruskal(gr))
    @test sum(kruskal(gr)) == sum(prim(gr))
    @test sum(kruskal(gr,kind=:min)) <= sum(kruskal(gr,kind=:max))
end

# every single function/method should have at least one test!

