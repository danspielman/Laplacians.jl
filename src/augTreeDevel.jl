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
    
