#=

This is a script that we use to run speed tests on Chimera graphs.
Run like julia st_chimera.jl 123456 10
Which will run for approximately 10 hours on graphs with 123456 vertices

=#

n = Base.parse(Int64,ARGS[1])
hours = Base.parse(Float64,ARGS[2])

fn = "st_chimera_onlyjl_$(n)h$(hours).jld"
println(fn)

using Laplacians
using MATLAB
using SparseArrays
using Statistics
using LinearAlgebra

lapdir = dirname(pathof(Laplacians))

include("$(lapdir)/../compare/matlabSafe.jl")
include("$(lapdir)/../compare/trilinosDrivers.jl")
include("$(lapdir)/../compare/compare_solvers_TL.jl")


ac_deg = function(a; verbose=false, args...)
    approxchol_lap(a; params=ApproxCholParams(:deg), verbose=verbose, args...)
end
ac_wdeg = function(a; verbose=false, args...)
    approxchol_lap(a; params=ApproxCholParams(:wdeg), verbose=verbose, args...)
end


test_ac = SolverTest(ac_deg, "ac")
test_acw = SolverTest(ac_wdeg, "ac_fast")
#test_amg = SolverTest(AMGLapSolver, "pyamg")
    
# removed chol because killing it can cause a crash
#test3 = SolverTest(chol_lap, "chol")
#test4 = SolverTest(cgLapSolver, "cg")
#tests = [test_ac test_acw test_amg]
#tests = [test_ac test_acw]

    tests = [test_acw]
    
# warm up the code
a = chimera(10000,1)
b = randn(10000)
b = b .- mean(b)
for solver in tests
    f = solver.solver(a,verbose=true)
    x = f(b)
end

dic = Dict()

using JLD2

i = 0
t0 = time()

while time() - t0 < 60*60*hours

    ti = time()
    global i += 1
    println("-----------")
    println(i)

    @time a = chimera(n,i)
    @time aw = randWeight(a)

    tn = "chimera($n,$i)"
    b = randn(n)
    b = b .- mean(b)
    b = b / norm(b)
    x = testVMatlabLap(tests, dic, a, b, verbose=true, tol=1e-8, testName=tn, test_icc=false, test_cmg=false, test_lamg=false, test_muelubelos=false )
    @save fn dic
    tn = "wtedChimera($n,$i)"
    x = testVMatlabLap(tests, dic, aw, b, verbose=true, tol=1e-8, testName=tn, test_icc=false, test_cmg=false, test_lamg=false, test_muelubelos=false )
    @save fn dic

    @printf("time this iter = %f",time() - ti)
    @printf("time all  iter = %f",time() - t0)

end
