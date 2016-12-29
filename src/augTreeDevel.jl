"""
    f = sddmWrapLap(lapSolver, sddm::AbstractArray; tol::Real=1e-6, maxits=Inf, maxtime=Inf, verbose=false, pcgIts=Int[], params...)
    f = sddmWrapLap(lapSolver)

Uses a `lapSolver` to solve systems of linear equations in sddm matrices.
"""
function sddmWrapLap(lapSolver, sddm::AbstractArray; tol::Real=1e-6, maxits=Inf, maxtime=Inf, verbose=false, pcgIts=Int[], params...)

    # Make a new adj matrix, a1, with an extra entry at the end.
    a, d = adj(sddm)
    a1 = extendMatrix(a,d)
    F = lapSolver(a1, tol=tol, maxits=maxits, maxtime=maxtime, verbose=verbose, pcgIts=pcgIts, params...)

    # make a function that solves the extended system, modulo the last entry
    tol_=tol
    maxits_=maxits
    maxtime_=maxtime
    verbose_=verbose
    pcgIts_=pcgIts

    f = function(b; tol=tol_, maxits=maxits_, maxtime=maxtime_, verbose=verbose_, pcgIts=pcgIts_)
        @time xaug = F([b; -sum(b)], tol=tol, maxits=maxits, maxtime=maxtime, verbose=verbose, pcgIts=pcgIts)
        xaug = xaug - xaug[end]
        return xaug[1:a.n]
    end
        
    return f
                                     
end

import Base.SparseArrays

macro cholmod_name(nm,typ) string("cholmod_", eval(typ) == ba.SuiteSparse_long ? "l_" : "", nm) end

function print_common{Tv<:ba.VTypes}(F::ba.Factor{Tv}, name::String)
    ba = Base.SparseArrays.CHOLMOD
    cm = ba.common()
        
    cm = ba.common()
    ba.set_print_level(cm, 3)
    ccall((@cholmod_name("print_common", ba.SuiteSparse_long),:libcholmod), 
    Cint, (Ptr{UInt8}, Ptr{UInt8}), name,  cm)
    nothing
end

"""
  nnzL, flops = ask_cholmod(mat)

Estimate the number of nonzeros in the cholfact factorization of mat, 
along with the number of flops needed to compute it.
Does this through a call to the analyze routine of cholmod.
Note that this is much faster than actually computing the factorization
"""
function ask_cholmod(sdd)
    ba = Base.SparseArrays.CHOLMOD
    cm = ba.common()

    anal = ba.analyze(ba.Sparse(lap(sdd)), cm);

    s_anal = unsafe_load(get(anal.p))

    n = Int(s_anal.n)
    
    nnzL = 0
    flops = 0
    
    for i in 1:n
        nzl = Int(unsafe_load(s_anal.ColCount,i))
        nnzL += nzl
        flops += nzl^2
    end

    return nnzL, flops
end
    
    
function print_cholmod(sdd)
    anal = ba.analyze(ba.Sparse(lap(sdd)), cm);

    print_common(anal, "F")

end


function augmentTree2{Tv,Ti}(tree::SparseMatrixCSC{Tv,Ti}, A::SparseMatrixCSC{Tv,Ti}, k::Ti)


    Aminus = A - tree
    
    ai,aj,av = findnz(triu(Aminus))
    m = length(ai)

    if m <= k
        return A
    end

    st = compStretches(tree, Aminus)
    _,_,sv = findnz(triu(st))

    
    r = -log(rand(m)) ./ sv
    ord = sortperm(r)

    edgeinds = ord[1:k]
    augi = ai[edgeinds]
    augj = aj[edgeinds]
    augv = av[edgeinds]
    

    n = size(A,1)
    aug = sparse(augi, augj, augv, n, n)
    aug = aug + aug'

    return tree + aug
    
end

function augmentTree3{Tv,Ti}(tree::SparseMatrixCSC{Tv,Ti}, A::SparseMatrixCSC{Tv,Ti}; nnzL_fac=2.0, flops_fac=100.0)

    Aminus = A - tree
    
    ai,aj,av = findnz(triu(Aminus))

    n = A.n
    m = length(ai)

    st = compStretches(tree, Aminus)
    _,_,sv = findnz(triu(st))

    
    r = -log(rand(m)) ./ sv
    ord = sortperm(r)

    #=
    @show n*nnzL_fac
    @show n*nnzL_fac/2
    @show (n+m)*flops_fac/4
    @show (n+m)*flops_fac
    =#
    
    nnzLTooBig(nnzL) = (nnzL-2*(n-1)) > n*nnzL_fac
    nnzLTooSmall(nnzL) = (nnzL-2*(n-1)) < n*nnzL_fac/2
    flopsTooBig(flops) = flops > (n+m)*flops_fac
    flopsTooSmall(flops) = flops < (n+m)*flops_fac/4
    
    k = Int(round(2*sqrt(n)))

    first = true
    direction = 0
    done = false

    while ~done
    
        edgeinds = ord[1:min(k,m)]
        augi = ai[edgeinds]
        augj = aj[edgeinds]
        augv = av[edgeinds]
        

        n = size(A,1)
        aug = sparse(augi, augj, augv, n, n)
        aug = aug + aug'

        augTree = tree+aug

        nnzL, flops = ask_cholmod(lap(augTree))

        @show k, nnzL, flops

        #=
        @show [n*nnzL_fac/2, nnzL-2*(n-1), n*nnzL_fac]
        @show [(n+m)*flops_fac/4, flops, (n+m)*flops_fac]
        =#
        
        if first
            first = false
            if ~(nnzLTooBig(nnzL) || nnzLTooSmall(nnzL) ||
                 flopsTooBig(flops) || flopsTooSmall(flops))
                done = true
            elseif nnzLTooBig(nnzL) || flopsTooBig(flops)
                k = div(k,2)
                direction = -1
            else
                k = k * 2
                direction = 1
                if k >= m
                    done = true
                end
                
            end
        else
            if direction == -1
                if nnzLTooBig(nnzL) || flopsTooBig(flops)
                    k = div(k,2)
                else
                    done = true
                end
                
            else # direction == 1

                if nnzLTooSmall(nnzL) && flopsTooSmall(flops)
                    k = k * 2
                    if k >= m
                        done = true
                    end
                else
                    @show nnzLTooSmall(nnzL)
                    @show flopsTooSmall(flops)
                    done = true
                end
            end
        end

        if done
            return augTree
        end
    end
        
    
end


function augTreeLapSolver3{Tv,Ti}(a::SparseMatrixCSC{Tv,Ti}; tol::Real=1e-6, maxits=Inf, maxtime=Inf, verbose=false, pcgIts=Int[], nnzL_fac=2.0, flops_fac=100.0 )

    t = akpw(a)
    at = augmentTree3(t,a; nnzL_fac=nnzL_fac, flops_fac=flops_fac)

    return pcgLapSolver(a, at, verbose=verbose, tol=tol, maxits=maxits, maxtime=maxtime, pcgIts=pcgIts)

end


#==========================================================
Speed Testing
===========================================================#
    

using DataFrames

global speedDF = DataFrame()


type SolverTest
    solver::Function
    name::String
end



function initDictCol!(dic, name, typ)
    if ~haskey(dic,name)
        dic[name] = Array(typ,0)
    end
end



"""

    function speedTestLapSolvers{Tv,Ti}(solvers, dic, a::SparseMatrixCSC{Tv,Ti}, b::Array{Tv,1}; tol::Real=1e-6, maxits=Inf, maxtime=Inf, verbose=false)

Runs many Laplacians solvers.  Puts the build and solve time results into solverDF.  Returns the answer given by KMP.

"""
function speedTestLapSolvers{Tv,Ti}(solvers, dic, a::SparseMatrixCSC{Tv,Ti}, b::Array{Tv,1}; tol::Real=1e-2, maxits=Inf, maxtime=Inf, verbose=false)


    b = b - mean(b)
    
    la = lap(a)

    it = Int[1]

    # init the dict, if it is the first time
    initDictCol!(dic, "nv", Int)
    initDictCol!(dic, "ne", Int)
    initDictCol!(dic, "hash_a", UInt64)
    
    solvecol(name) = "$(name)_solve"
    buildcol(name) = "$(name)_build"
    totcol(name) = "$(name)_tot"
    itscol(name) = "$(name)_its"
    errcol(name) = "$(name)_err"
    
    for solver in solvers
        name = solver.name
        initDictCol!(dic, solvecol(name), Float64)
        initDictCol!(dic, buildcol(name), Float64)
        initDictCol!(dic, totcol(name), Float64)
        initDictCol!(dic, itscol(name), Float64)
        initDictCol!(dic, errcol(name), Float64)
    end
    
    nv = size(a,1)
    ne = nnz(a)
    hash_a = hash(a)

    push!(dic["nv"],nv)
    push!(dic["ne"],ne)
    push!(dic["hash_a"],hash_a)
    
    for solverTest in solvers
        tic()
        f = solverTest.solver(a, tol=tol, maxtime=maxtime, maxits=maxits, verbose=verbose)
        build_time = toq()

        tic()
        x = f(b, pcgIts = it)
        solve_time = toq()

        err = norm(la * x - b) / norm(b)

        name = solverTest.name
        push!(dic[solvecol(name)],solve_time)
        push!(dic[buildcol(name)],build_time)
        push!(dic[totcol(name)],solve_time + build_time)
        push!(dic[itscol(name)],it[1])
        push!(dic[errcol(name)],err)

    end

    return x

end



