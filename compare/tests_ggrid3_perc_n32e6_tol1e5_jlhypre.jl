#=

This is a script that we use to run speed tests on two-dim grids.
Run like julia st_grid3.jl 1
Which will run for approximately 1 hour, increasing the sizes of the graphs.

=#

hours = Base.parse(Float64,ARGS[1])

fn = "$(PROGRAM_FILE).h$(hours).jld"
println(fn)

using Laplacians
using MATLAB
using SparseArrays
using Statistics
using LinearAlgebra

lapdir = dirname(pathof(Laplacians))

include("$(lapdir)/../compare/matlabSafe.jl")
include("$(lapdir)/../compare/hypreDrivers.jl")
include("$(lapdir)/../extern/hypre/hypreExport.jl")
include("$(lapdir)/../compare/compare_solvers_TL.jl")

ac_deg = function(a; verbose=false, args...)
    approxchol_sddm(a; params=ApproxCholParams(:deg), verbose=verbose, args...)
end
# ac_wdeg = function(a; verbose=false, args...)
#     approxchol_sddm(a; params=ApproxCholParams(:wdeg), verbose=verbose, args...)
# end


test_ac = SolverTest(ac_deg, "ac")
#test_acw = SolverTest(ac_wdeg, "ac_fast")
#test_amg = SolverTest(AMGLapSolver, "pyamg")
#test_chol = SolverTest(chol_sddm, "chol") 
#tests = [test_ac test_acw test_amg test_chol]
#tests = [test_ac test_acw test_chol]
#tests = [test_ac test_acw]
tests = [test_ac]

using JLD2

# warm up the test code
dicWarmup = Dict()
n = 10000
println("----- warm up starting ------")
println(n)

s = round(Int,n^(1/3))
s1 = s
s2 = s
s3 = s
w = 1e5
seed = 1234
tn = "ggrid3_perc($(s1),$(s2),$(s3),$(w),$(seed))"
@time A = ggrid3_perc(s1,s2,s3,w,seed);
@time L = lap(A)
@time int = getInterior3(s1,s2,s3)
@time M = L[int,int];
ni = length(int)
b = randn(ni);
b = b / norm(b)
x = testSddm(tests, dicWarmup, M, b; verbose=false, tol=1e-5, testName=tn, test_icc=false, test_cmg=false, test_lamg=false)
@save fn dicWarmup    

println("----- warm up complete ------")

# actual tests

using JLD2
dic = Dict()

n = round(Int,32e6)
w = 1e5

iter = 0

t0 = time()

seed = 1234+32
# I think n=32e6 and seed = 1234+33 was where we first had trouble? 
# note that seed is incremented before use

while time() - t0 < 60*60*hours

    println("-----------")
    println("n = $(n)")

    s = round(Int,n^(1/3))
    s1 = s
    s2 = s
    s3 = s
    global seed += 1
    tn = "ggrid3_perc($(s1),$(s2),$(s3),$(w),$(seed))"
    @time A = ggrid3_perc(s1,s2,s3,w,seed);
    @time L = lap(A)
    @time int = getInterior3(s1,s2,s3)
    @time M = L[int,int];
    ni = length(int)
    b = randn(ni);
    b = b / norm(b)
    x = testSddm(tests, dic, M, b; verbose=true, tol=1e-5, testName=tn, test_icc=false, test_cmg=false, test_lamg=false)
    @save fn dic

    global iter += 1
end
