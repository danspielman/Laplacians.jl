n = parse(Int,ARGS[1])
hours = parse(Float64,ARGS[2])

t_start = time()

jld = "mcfp_chim_$(n)h$(hours).jld"
println(jld)

using Laplacians

lapdir = Pkg.dir("Laplacians")
include("$(lapdir)/src/flowUtils.jl")
include("$(lapdir)/src/min_cost_flow.jl")
include("$(lapdir)/src/flow_extern.jl")


# warm up run
k = round(Int,n/20)
mcfp = gen_mcfp(wtedChimera(n,1), frac_of_max=0.9, k=k, choice = :spect);
fn = "tmp.txt"
writeDimacsMCF(fn,mcfp)

cs2_time, _ = call_cs2(fn)
lemon_time, _ = call_lemon(fn)
t1 = time()
sol = min_cost_flow(Float(mcfp), lapSolver = h->approxCholSddm(h,tol=1e-8), tol = 1e-5)
ipm_time = time()-t1


dic = Dict()
dic["cs2"] = []
dic["lemon"] = []
dic["ipm"] = []


using JLD

i = 0

while time() - t_start < 60*60*hours

    i += 1
    println("----------------------")
    println(i)

    mcfp = gen_mcfp(wtedChimera(n,i), frac_of_max=0.9, k=k, choice = :spect);
    fn = "tmp.txt"
    writeDimacsMCF(fn,mcfp)

    @show now()
    cs2_time, _ = call_cs2(fn)
    @show now()
    lemon_time, _ = call_lemon(fn)
    @show now()
    t1 = time()
    sol = min_cost_flow(Float(mcfp), lapSolver = h->approxCholSddm(h,tol=1e-8), tol = 1e-5)
    ipm_time = time()-t1

    rm(fn)

    push!(dic["lemon"], lemon_time)
    push!(dic["cs2"], cs2_time)
    push!(dic["ipm"], ipm_time)

    save(jld,"dic",dic)

end
