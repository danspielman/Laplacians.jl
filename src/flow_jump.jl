using JuMP
using Clp

function stFlow_jump(stfp::stFlowProblem)
    edge_list = stfp.edge_list
    m = size(edge_list,1)
    n = maximum(edge_list)

#    internal = ones(Bool,n)
#    internal[stfp.s] = false
#    internal[stfp.t] = false
    
    B = sparse(collect(1:m), edge_list[:,1], 1.0, m, n) -
      sparse(collect(1:m), edge_list[:,2], 1.0, m, n)
    B[:,stfp.s] = 0
    B[:,stfp.t] = 0

    mod = Model(solver=ClpSolver())
    @variable(mod, x[1:m] >= 0)
    @constraint(mod, x .<= stfp.capacities)
    @constraint(mod, B'*x .== 0.0)
    obj = zeros(m)
    idx = find(edge_list[:,1].==stfp.s)
    obj[idx] = 1

    @objective(mod, Max, obj'*x)
    @time status = solve(mod)

    f = getvalue(x)
    return f
    
end



function MCFjump(mcfp::MCFproblem)
    edge_list = mcfp.edge_list
    m = size(edge_list,1)
    n = maximum(edge_list)
    B = sparse(collect(1:m), edge_list[:,1], 1.0, m, n) -
      sparse(collect(1:m), edge_list[:,2], 1.0, m, n)
 

    mod = Model(solver=ClpSolver())
    @variable(mod, x[1:m] >= 0)
    @constraint(mod, x .<= mcfp.capacities)
    @constraint(mod, B'*x .== mcfp.demands)
    @objective(mod, Min, mcfp.costs'*x)
    @time status = solve(mod)

    f = getvalue(x)
    return f
    
end

