

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
  mcfp = readDimacsMCF(filename)

Read a minimum cost flow problem in the <a href="http://lpsolve.sourceforge.net/5.5/DIMACS_mcf.htm">Dimacs format</a>.
Our code only handles 0 lower bounds on flows.
"""
function readDimacsMCF(filename)

    fi = open(filename)
    lines = readlines(fi)
    close(fi)

    edge_list = []
    dems = []
    costs = []
    caps = []
    
    edge_list_ptr = 0
    for line in lines

        if line[1] == 'p'
            line_parts = split(line)
            n = parse(Int,line_parts[3])
            m = parse(Int,line_parts[4])
            edge_list = zeros(Int,m,2)
            caps = zeros(Float64,m)
            costs = zeros(Float64,m)
            dems = zeros(Float64,n)
        end

        if line[1] == 'n'
            line_parts = split(line)
            u = parse(Int,line_parts[2])
            d = parse(Float64,line_parts[3])
            dems[u] = d 
        end
        
        if line[1] == 'a'
            line_parts = split(line)
            src = parse(Int,line_parts[2])
            dst = parse(Int,line_parts[3])
            low = parse(Int,line_parts[4])
            cap = parse(Float64,line_parts[5])
            cost = parse(Float64,line_parts[6])
            

            if (low > 0) 
                error("we only handle zero lower bounds on flows")
            end
            edge_list_ptr += 1
            edge_list[edge_list_ptr,:] = [src dst]
            caps[edge_list_ptr] = cap
            costs[edge_list_ptr] = cost
        end
    end

    mcfp = MCFproblem(edge_list, caps, costs, dems)

    return mcfp
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




function makeSTFlowProblem(a; k=round(Int,sqrt(size(a,1))))

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
    order = Array(Ti,n)

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

