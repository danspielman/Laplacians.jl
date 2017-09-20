#=

NOTE: calls to JuMP are now in flow_jump.jl
Function related to Dimacs files and other solvers are now in flow_extern.jl

Include with the lines
lapdir = Pkg.dir("Laplacians")
include("$(lapdir)/src/flowUtils.jl")

Contains:
function reportMCFresults(mcfp, flow)
function check_flow(stfp::stFlowProblem, f; tol=1e-3)
function quant_vec(v, maxint)
prob = makeSTFlowProblem(a; k=round(Int,sqrt(size(a,1))), choice=:farpair)
function makeSTFlowProblemOld(a; k=round(Int,sqrt(size(a,1))))
function Float(mcfp::MCFproblem)
function rand_mcfp(a;  maxcost=1000, maxcap=1000)
function gen_mcfp(a;
mcfp = gen_mcfp(a; maxcost=, maxcap=, frac_of_max=,k=,choice=)
function reachFrom{Tv,Ti}(dirg::SparseMatrixCSC{Tv,Ti}, v0::Ti)
function cutCapacity{Tv,Ti}(dirg::SparseMatrixCSC{Tv,Ti}, comp::Array{Ti,1})


=#

"""
type MCFproblem{Tv,Ti}
    edge_vertex_mat::SparseMatrixCSC{Tv,Ti}
    costs::Array{Tv,1}
    capacities::Array{Tv,1}
    demands::Array{Tv,1}
end

  Each row of the edge_list has the origin followed by the destination of the edge.
  The positive demands are the on the nodes that flow comes from, and it goes to negative.
  So, if there is just a source and sink, then the demand at the source should be positive.
"""
type MCFproblem{Tv,Ti}
    edge_list::Array{Ti,2}
    capacities::Array{Tv,1}
    costs::Array{Tv,1}
    demands::Array{Tv,1}
end


"""
  A Max Flow Problem is specified by an edge list (from to) for each,
  capacities on edges,
  and a demands on vertices that sum to 0
"""
type MaxFlowProblem{Tv,Ti}
    edge_list::Array{Ti,2}
    capacities::Array{Tv,1}
    demands::Array{Tv,1}
end

"""
  An st Flow Problem is specified by an edge list (from to) for each,
  capacities on edges,
  s and t
"""
type stFlowProblem{Tv,Ti}
    edge_list::Array{Ti,2}
    capacities::Array{Tv,1}
    s::Int
    t::Int
end



"""
    cost = reportMCFresults(mcfp, flow)

Given an min cost flow problem, and a solution flow, report how well it does
and return the cost.
"""
function reportMCFresults(mcfp, flow)
    cost = mcfp.costs'*flow
    println("Cost: ", cost)
    println("Min flow: ", minimum(flow))
    println("Min slack: ", minimum(mcfp.capacities-flow))

    edge_list = mcfp.edge_list
    m = size(edge_list,1)
    n = maximum(edge_list)
    B = sparse(collect(1:m), edge_list[:,1], 1.0, m, n) -
    sparse(collect(1:m), edge_list[:,2], 1.0, m, n)

    println("Error on demands: ", sum(abs.(B'*flow- mcfp.demands)))

    return cost
end

"""
    flowval = check_flow(stfp::stFlowProblem, f; tol=1e-3)

    Given an st flow problem and a flow f, checks how good the flow is,
    and returns its value.
    It constructs a dual cut, using all edges flowing at least 1-tol times their capacity.
"""
function check_flow(stfp::stFlowProblem, f; tol=1e-3)
    n = maximum(stfp.edge_list)
    m = size(stfp.edge_list,1)

    obj = zeros(m)
    idx = find(stfp.edge_list[:,1].==stfp.s)
    obj[idx] = 1

    value = obj'*f
    println("Value: ", value)

    unsaturated = f .< (1-tol)*stfp.capacities

    usei = stfp.edge_list[unsaturated,1]
    usej = stfp.edge_list[unsaturated,2]
    adir = sparse(usej,usei,1,n,n) # note backwards convention
    comp = reachFrom(adir,stfp.s)

    a = sparse(stfp.edge_list[:,1],stfp.edge_list[:,2],stfp.capacities,n,n)
    cut_val = cutCapacity(a,comp)

    println("The cut value is ", cut_val)

    scut = sum(stfp.capacities[stfp.edge_list[:,1] .== stfp.s])
    tcut = sum(stfp.capacities[stfp.edge_list[:,2] .== stfp.t])
    println("cut {s} has value $(scut) and cut {t} has value $(tcut)")

    if any(f.<0)
        println("most negative flow: ", minimum(f[f.<0]))
    end
    over = f-stfp.capacities
    if any(over.>0)
        println("most over capacity: ", maximum(over[over.>0]))
    end

    B = sparse(collect(1:m), stfp.edge_list[:,1], 1.0, m, n) -
      sparse(collect(1:m), stfp.edge_list[:,2], 1.0, m, n)
    B[:,stfp.s] = 0
    B[:,stfp.t] = 0
    outflow = f'*B
    println("flow-in = flow-out failure: ",maximum(abs.(outflow)))

    return value

end



"""
    v = quant_vec(x, maxint)

Quantize the vector x to integers of absolute value at most maxint.
"""
function quant_vec(v, maxint)
    maxv = maximum(abs.(v))
    if maxv == 0
        return round.(Int,v)
    elseif (maxv < maxint) && (maximum(v-round.(Int,v)) == 0)
        return round.(Int,v)
    else
        f(x) = round(Int,sign(x))*ceil(Int, maxint * abs(x) / maxv)
        return f.(v)
    end
end



"""
    prob = makeSTFlowProblem(a; k=round(Int,sqrt(size(a,1))), choice=:farpair)

Create an instance of an stFlowProblem with one source and once sink
from a graph given by its Adjacency matrix `a` where the weights are capacities.
It creates two new vertices s and t.
These are each connected to `k` vertices in the graph.
These k vertices are chosen according to `choice`.

* :farpair means choose one vertex and the k closest to it.  Then choose the vertex farthese from that, and the k closest to it.
* :spect means approximate a fiedler vector, use it to order, and choose the extreme k.
* :rand choose those k at random

"""
function makeSTFlowProblem(a; k=round(Int,sqrt(size(a,1))), choice=:farpair)

    n = size(a,1)

    if choice == :farpair

        s0 = rand(1:n)
        ds0 = shortestPaths(a,s0)[1]
        sk = sortperm(ds0)[1:k]

        t0 = findmax(ds0)[2]
        dt0 = shortestPaths(a,t0)[1]
        tk = sortperm(dt0)[1:k]


    elseif choice == :spect
        f = approxCholLap(a)
        op = SqLinOp(true,1.0,size(a,1),f)
        v = vec(eigs(op;nev=1,which=:LM)[2])
        p = sortperm(v)
        sk = p[1:k]
        tk = p[(n-k+1):n]

    else  # choice == :rand
        p = randperm(n)
        sk = p[1:k]
        tk = p[(k+1):(2*k)]
    end

    deg = sum(a,2)

    s = n+1
    t = n+2
    anew = spzeros(n+2,n+2)

    anew[1:n,1:n] = a
    # anew[sk,s] = deg[sk]
    anew[s,sk] = deg[sk]
    anew[tk,t] = deg[tk]
    # anew[t,tk] = deg[tk]

    (ai,aj,av) = findnz(anew)
    return stFlowProblem([ai aj],av,s,t)

end

function makeSTFlowProblemOld(a; k=round(Int,sqrt(size(a,1))))

    n = size(a,1)
    s0 = rand(1:n)
    ds0 = shortestPaths(a,s0)[1]
    sk = sortperm(ds0)[1:k]

    t0 = findmax(ds0)[2]
    dt0 = shortestPaths(a,t0)[1]
    tk = sortperm(dt0)[1:k]

    deg = sum(a,2)

    s = n+1
    t = n+2
    anew = spzeros(n+2,n+2)

    anew[1:n,1:n] = a
    anew[sk,s] = deg[sk]
    anew[s,sk] = deg[sk]
    anew[tk,t] = deg[tk]
    anew[t,tk] = deg[tk]

#=
    demand = min(sum(deg[tk]), sum(deg[sk]))

    BBt = -dirEdgeVertexMat(anew)'
    u = copy(anew.nzval)
    b = zeros(n+2)
    b[s]=demand
    b[t]=-demand
=#

    return anew, s, t
#    return BBt, u, b
end


"""
    mcfp = stf_to_mcf(stfp::stFlowProblem)

Recast an st flow problem as a min cost flow problem.
It adds one edge at the end, so to recover the st flow one needs to remove it
"""
function stf_to_mcf(stfp::stFlowProblem)

    n = maximum(stfp.edge_list)
    m = size(stfp.edge_list,1)

    edges = zeros(Int, m+1, 2)
    for i in 1:m
        edges[i,:] = stfp.edge_list[i,:]
    end
    edges[m+1,:] = [stfp.t stfp.s]

    costs = zeros(m+1)
    costs[m+1] = -1

    dems = zeros(n)
    caps = copy(stfp.capacities)

    push!(caps,sum(caps[find(edges[:,1].==stfp.s)]))

    mcfp = MCFproblem(edges, caps, costs, dems)
end


function Float(mcfp::MCFproblem)
    return MCFproblem(mcfp.edge_list, Float64.(mcfp.capacities), Float64.(mcfp.costs), Float64.(mcfp.demands))
end


"""
    function rand_mcfp(a; maxdem=1000, maxcost=1000, maxcap=1000)

Make a random Min Cost Flow problem from the graph `a`
by setting all demands to zero, the costs to random positive and negative numbers,
and the capacities to random integers between 0 and maxcap
"""
function rand_mcfp(a;  maxcost=1000, maxcap=1000)

    n = size(a,1)

    ai, aj, av = findnz(a)
    m = length(ai)
    caps = quant_vec(av, maxcap)

    dems = zeros(Int,n)

    costs = quant_vec( (rand(m)-0.5), maxcost)

    mcfp = MCFproblem([ai aj], caps, costs, dems)

end

"""
    mcfp = gen_mcfp(a; maxcost=, maxcap=, frac_of_max=,k=,choice=)

Generate an mcfp problem from the graph a.
First generate an st-flow problem, by adding vertices s and t.
Find the max flow in this graph (approximately),
and then ask to flow frac_of_max of that.
Generate the costs at random, and quantitze the capacities.
k and choice are inputs to the st-flow problem generator.
"""
function gen_mcfp(a;
                  maxcost=1000,
                  maxcap = 1000,
                  frac_of_max = 0.9,
                  k=round(Int,sqrt(size(a,1))),
                  choice=:spect)

    stfp = makeSTFlowProblem(a; k = k, choice = choice )

    stfp.capacities = quant_vec(stfp.capacities, maxcap)

    as_mcf = stf_to_mcf(stfp)
    sol = min_cost_flow(Float(as_mcf), lapSolver = approxCholSddm, tol = 1e-3)
    flow = sol[1][1:(end-1)]
    val = round(Int,check_flow(stfp, flow) * frac_of_max)

    m = size(stfp.edge_list,1)
    n = size(a,1)+2

    dems = zeros(n)
    dems[stfp.s] = val
    dems[stfp.t] = -val

    costs = Float64.(quant_vec(rand(m), maxcost))

    mcfp = MCFproblem(stfp.edge_list, stfp.capacities, costs, dems)

    return mcfp
end


function stFlow(dirg, s0, t0; lapSolver = cholLap, tol=1e-6)

    U = sum(dirg[:,s0])*2

    Bt = -dirEdgeVertexMat(dirg)'

    @show size(Bt)

    n,m = size(Bt)

    b = zeros(n)
    b[s0] = U
    b[t0] = -U

    cap = copy(dirg.nzval)
    # cap = ones(m)+0.0
    (x,y,s) = max_flow_IPM(Bt,b,cap,lapSolver = lapSolver,tol=tol)

    f = x[1:m]
    t = x[end]

    println()

    println("The error in satisfying the demand is ", norm(Bt*f-t*b))
    println("The flow value is ", t*U)

    (ai,aj) = findnz(dirg)
    use = f .< 0.95*cap

    usei = ai[use]
    usej = aj[use]
    adir = sparse(usei,usej,1,n,n)
    comp = reachFrom(adir,s0)

    cut_val = cutCapacity(dirg,comp)

    println("The cut value is ", cut_val)


    flow_val = x[end]*U

    return flow_val, cut_val, f, comp

end



"""
  edge dirg[i,j] = 1 is an edge from j to i
  Warning: this might look backwards
"""
function reachFrom{Tv,Ti}(dirg::SparseMatrixCSC{Tv,Ti}, v0::Ti)
    n = dirg.n

    comp = zeros(Ti,n)
    order = Array{Ti}(n)

    # note that all of this casting is unnecessary.
    # but, some of it speeds up the code
    # I have not figured out the minimal necessary
    c::Ti = 0

    colptr::Array{Ti,1} = dirg.colptr
    rowval::Array{Ti,1} = dirg.rowval


    comp[v0] = 1

    ptr::Ti = 1
    orderLen::Ti = 2
    order[ptr] = v0

    while ptr < orderLen
        curNode = order[ptr]

        for ind in colptr[curNode]:(colptr[curNode+1]-1)
            nbr = rowval[ind]
            if comp[nbr] == 0
                comp[nbr] = 1
                order[orderLen] = nbr
                orderLen += 1
            end # if
        end # for
        ptr += 1
    end # while

  return comp
end # function

"""
  edge dirg[i,j] = 1 is an edge from j to i
  Warning: this might look backwards
"""
function cutCapacity{Tv,Ti}(dirg::SparseMatrixCSC{Tv,Ti}, comp::Array{Ti,1})

    cap = 0

    colptr::Array{Ti,1} = dirg.colptr
    rowval::Array{Ti,1} = dirg.rowval

    for j in 1:dirg.n
        if comp[j] == 1
            for i in colptr[j]:(colptr[j+1]-1)
                if comp[rowval[i]] == 0
                    cap += dirg.nzval[i]
                end
            end
        end
    end
    return cap

end
