#=

This is a script that we use to run speed tests on two-dim grids.
Run like julia st_grid3.jl 1
Which will run for approximately 1 hour, increasing the sizes of the graphs.

=#

hours = Base.parse(Float64,ARGS[1])

fn = "tests_ggrid3_h$(hours).jld"
println(fn)

using Laplacians
using MATLAB
using SparseArrays
using Statistics
using LinearAlgebra

lapdir = dirname(pathof(Laplacians))

include("$(lapdir)/../compare/matlabSafe.jl")
include("$(lapdir)/../compare/hypreDrivers.jl")
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
test_chol = SolverTest(chol_sddm, "chol")
#tests = [test_ac test_acw test_amg test_chol]
tests = [test_ac test_acw test_chol]
#tests = [test_ac test_acw]

# warm up the code
n = 10000
s = round(Int,n^(1/3))
s1 = s
s2 = s
s3 = s
@time A = ggrid3(s1,s2,s3);
@time L = lap(A)
@time int = getInterior3(s1,s2,s3)
@time M = L[int,int];
tn = "ggrid3($(s))"
ni = length(int)
b = randn(ni);
b = b / norm(b)
for solver in tests
    f = solver.solver(M,verbose=true)
    x = f(b)
end

using JLD2
dic = Dict()



n = 10000
t0 = time()

while time() - t0 < 60*60*hours

    println("-----------")
    println(n)

    s = round(Int,n^(1/3))
    s1 = s
    s2 = s
    s3 = s
    @time A = ggrid3(s1,s2,s3);
    @time L = lap(A)
    @time int = getInterior3(s1,s2,s3)
    @time M = L[int,int];
    tn = "ggrid3($(s))"
    ni = length(int)
    b = randn(ni);
    b = b / norm(b)
    x = testSddm(tests, dic, M, b; verbose=false, tol=1e-8, testName=tn, test_icc=true, test_cmg=false, test_lamg=false)
    @save fn dic
    global n = 2*n

end
