#=

Code for using external flow problem generators and solvers,
calling them from the shell

Include with the lines
lapdir = Pkg.dir("Laplacians")
include("$(lapdir)/src/flow_extern.jl")

Provides
  mcfp = readDimacsMCF(filename)
  function writeDimacsMCF(filename, mcfp)
  function goto8_gen(n, seed; maxcost = 1000, maxcap = 1000)
  function call_cs2(fn)
  function call_lemon(fn)


These need compiled programs:
goto
cs2
lemon (which is a link to dimacs_solver)
=#


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
   writeDimacsMCF(filename, mcfp)

Write a minimum cost flow problem in the <a href="http://lpsolve.sourceforge.net/5.5/DIMACS_mcf.htm">Dimacs format</a>.
"""
function writeDimacsMCF(filename, mcfp)

    fi = open(filename, "w")

    edge_list = mcfp.edge_list
    m = size(edge_list,1)
    n = maximum(edge_list)

    # the problem line
    write(fi, "p min $(n) $(m)\n")


    dems = round.(Int32, quant_vec(mcfp.demands, round(Int(2^31-1))))
    for i in 1:n
        write(fi, "n $(i) $(dems[i])\n")
    end

    caps = round.(Int32, quant_vec(mcfp.capacities, 2^31-1))
    costs = round.(Int32, quant_vec(mcfp.costs, 2^31-1))

    for i in 1:m
        write(fi, "a $(edge_list[i,1]) $(edge_list[i,2]) 0 $(caps[i]) $(costs[i])\n")
    end

    close(fi)
end

"""
    goto8_gen(n, seed; maxcost = 1000, maxcap = 1000)

Call Goldberg's GOTO generator with m = 8*n, as per ...,
Set seed_base to 0 if you want to set the seed to goto exactly.
"""
function goto8_gen(n, seed; maxcost = 1000, maxcap = 1000, seed_base = 24369)
    goto_gen(n, 8*n, seed, maxcost=maxcost, maxcap=maxcap, seed_base=seed_base)
end


"""
    goto_gen(n, m, seed; maxcost = 1000, maxcap = 1000)

Call Goldberg's GOTO generator 
Set seed_base to 0 if you want to set the seed to goto exactly.
"""
function goto_gen(n, m, seed; maxcost = 1000, maxcap = 1000, seed_base = 24369)

    @assert n >= 15
    @assert m >= 6*n
    @assert m <= n^(5/3)
    @assert maxcap >= 8
    @assert maxcost >= 8

    fh = open("params.txt","w")
    write(fh,"$(n) $(m) $(maxcap) $(maxcost) $(seed+seed_base)")
    close(fh)

    fn = "goto8_$(n)_$(seed).min"
    run(pipeline( `goto`, stdin="params.txt", stdout=fn))

    return fn
end


"""
    runtime, cost = call_cs2(fn)

Call Goldberg's cs2 code on file fn, returning time and cost
"""
function call_cs2(fn)
    run(pipeline( `cs2`, stdin=fn, stdout="out.txt"))
    fh = open("out.txt")

    lines = readlines(fh)
    line = ""
    for line in lines
        if length(line) > 6 && line[1:6] == "c time"
            break
        end
    end

    parts = split(line)
    return parse(Float64,parts[3]), parse(Float64, parts[5])
end

"""
    runtime, cost = lemon(fn)

Lemon's MCF code on file fn, returning time and cost.
"""
function call_lemon(fn)
    run(pipeline( `lemon $(fn)`,  stderr="out.txt"))
    fh = open("out.txt")

    lines = readlines(fh)
    t = parse(Float64,split(lines[6])[end][1:(end-1)])
    c = parse(Int,split(lines[9])[end])

    return t, c
end
