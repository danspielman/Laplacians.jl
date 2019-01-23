#=

This is a script that we use to run speed tests on Chimera graphs.
Run like julia st_chimera.jl 123456 10
Which will run for approximately 10 hours on graphs with 123456 vertices

=#

# tprepre = time()

# __precompile()__

# @printf("precompile time = %.1f\n",tpre-tprepre)


tpre = time()

n = Base.parse(Int64,ARGS[1])
hours = Base.parse(Float64,ARGS[2])

fn = "st_chimera_onlymb_$(n)h$(hours).jld"
println(fn)

using Laplacians
using MATLAB
using SparseArrays
using Statistics
using LinearAlgebra
using Printf

#diagnostics
# using InteractiveUtils
# using Profile

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
tests = [test_ac test_acw]
#    tests = [test_acw]

run_icc = true
run_cmg = true
run_lamg = true
run_muelubelos = true
    
# warm up the solver code
a = chimera(1000,1)
b = randn(1000)
b = b .- mean(b)
for solver in tests
    f = solver.solver(a,verbose=true)
    x = f(b)
end

# warm up the test code
using JLD2
tn = "chimera(1000,1)"
dicWarmup = Dict() 
    @show @elapsed x = testVMatlabLap(tests, dicWarmup, a, b, verbose=true, tol=1e-8, testName=tn, test_icc=run_icc, test_cmg=run_cmg, test_lamg=run_lamg, test_muelubelos=run_muelubelos )
    @show @elapsed @save "warmupdummy.jld" dicWarmup
    
dic = Dict()




i = 0
t0 = time()
@printf("test load&warmup time = %.1f\n",t0-tpre)

while time() - t0 < 60*60*hours

    ti = time()
    global i += 1
    println("-----------")
    println(i)

    @time a = chimera(n,i)
    @show Base.summarysize(a)
    @time aw = rand_weight(a)
    @show Base.summarysize(aw)

    tn = "chimera($n,$i)"
    @time b = randn(n)
    @time b = b .- mean(b)
    @time b = b / norm(b)
    x = testVMatlabLap(tests, dic, a, b, verbose=true, tol=1e-8, testName=tn, test_icc=run_icc, test_cmg=run_cmg, test_lamg=run_lamg, test_muelubelos=run_muelubelos )
    @save fn dic
    tn = "wtedChimera($n,$i)"
    x = testVMatlabLap(tests, dic, aw, b, verbose=true, tol=1e-8, testName=tn, test_icc=run_icc, test_cmg=run_cmg, test_lamg=run_lamg, test_muelubelos=run_muelubelos )
    @save fn dic

    @printf("time (sec) this iter = %.1f\n",time() - ti)
    @printf("time (sec) all  iter = %.1f\n",time() - t0)

end

# li, lidict = Profile.retrieve()
# using JLD
# @save "/tmp/foo.jlprof" li lidict
