#=

This is a script that we use to run speed tests on two-dim grids.
Run like julia scriptname.jl 2000 1
Which will run for approximately 1 hour, 
on chimeras with n = 2000
=#

tpre = time()

n = 0

try
    global n
    n = convert(Int64,Base.parse(Float64,ARGS[1]))
catch e
    global n
    n = Base.parse(Int64,ARGS[1])
end

hours = Base.parse(Float64,ARGS[2])

fn = "$(PROGRAM_FILE).n$(n).h$(hours).jld"
println(fn)

using Laplacians
using MATLAB
using SparseArrays
using Statistics
using LinearAlgebra
using Printf

lapdir = dirname(pathof(Laplacians))

include("$(lapdir)/../compare/matlabSafe.jl")
include("$(lapdir)/../compare/hypreDrivers.jl")
include("$(lapdir)/../extern/hypre/hypreExport.jl")
include("$(lapdir)/../compare/compare_solvers_TL.jl")

ac_deg = function(a; verbose=false, args...)
    approxchol_sddm(a; params=ApproxCholParams(:deg), verbose=verbose, args...)
end
ac_wdeg = function(a; verbose=false, args...)
    approxchol_sddm(a; params=ApproxCholParams(:wdeg), verbose=verbose, args...)
end


test_ac = SolverTest(ac_deg, "ac")
test_acw = SolverTest(ac_wdeg, "ac_fast")
#test_amg = SolverTest(AMGLapSolver, "pyamg")
#test_chol = SolverTest(chol_sddm, "chol")
#tests = [test_ac test_acw test_amg test_chol]
#tests = [test_ac test_acw test_chol]
tests = [test_ac test_acw]

run_hypre = true
run_icc = false
run_cmg = false
run_lamg = false
run_muelubelos = false

using JLD2

# warm up the test code
println("----- warm up starting ------")
dicWarmup = Dict()
iWarmup = 0
nWarmup = 1000
println("i = $(iWarmup)")
println("n = $(nWarmup)")

@time a = chimera(nWarmup,iWarmup)
@show Base.summarysize(a)
L = lap(a)
@time aw = rand_weight(a)
@show Base.summarysize(aw)
Lw = lap(aw)

@time int = setdiff(1:nWarmup,1:ceil(nWarmup^(1/3)):nWarmup)
ni = length(int)
b = randn(ni);
@time b = b / norm(b);
@time M = L[int,int];
@time Mw = Lw[int,int];

tn = "chimera($nWarmup,$iWarmup)"
x = testVMatlabLap(tests, dicWarmup, M, b, verbose=true, tol=1e-8, testName=tn, test_icc=run_icc, test_cmg=run_cmg, test_lamg=run_lamg, test_hypre=run_hypre )
@save fn dicWarmup
tn = "wtedChimera($nWarmup,$iWarmup)"
x = testVMatlabLap(tests, dicWarmup, Mw, b, verbose=true, tol=1e-8, testName=tn, test_icc=run_icc, test_cmg=run_cmg, test_lamg=run_lamg, test_hypre=run_hypre )
@save fn dicWarmup


println("----- warm up complete ------")

# actual tests

using JLD2
dic = Dict()

i = 0
t0 = time()
@printf("test load&warmup time = %.1f\n",t0-tpre)

while time() - t0 < 60*60*hours

    ti = time()
    global i += 1
    println("-----------")
    println("i = $(i)")
    println("n = $(n)")

    @time a = chimera(n,i)
    @show Base.summarysize(a)
    L = lap(a)
    @time aw = rand_weight(a)
    @show Base.summarysize(aw)
    Lw = lap(aw)

    @time int = setdiff(1:n,1:ceil(n^(1/3)):n)
    ni = length(int)
    b = randn(ni);
    @time b = b / norm(b);
    @time M = L[int,int];
    @time Mw = Lw[int,int];

    tn = "chimera($n,$i)"
    x = testSddm(tests, dic, M, b, verbose=true, tol=1e-8, testName=tn, test_icc=run_icc, test_cmg=run_cmg, test_lamg=run_lamg, test_hypre=run_hypre )
    @save fn dic
    tn = "wtedChimera($n,$i)"
    x = testSddm(tests, dic, Mw, b, verbose=true, tol=1e-8, testName=tn, test_icc=run_icc, test_cmg=run_cmg, test_lamg=run_lamg, test_hypre=run_hypre )
    @save fn dic

    @printf("time (sec) this iter = %.1f\n",time() - ti)
    @printf("time (sec) all  iter = %.1f\n",time() - t0)


end
