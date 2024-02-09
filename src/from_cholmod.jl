#=
  Routines for finding out what cholmod would do when Cholesky factorizing a matrix

  These are calling cholmod's analyze routinem which is overkill.
  They could be made faster with more direct access to Cholmod.

  By Daniel Spielman.
=#

function get_suitesparse_handle()
    local cholmod = SuiteSparse.CHOLMOD
    local common =
        (VERSION.major >= 1 && VERSION.minor >= 9) ? cholmod.getcommon() :
        cholmod.COMMONS[Threads.threadid()]
    return (cholmod, common)
end


"""
  nnzL, flops = ask_cholmod(mat)

Estimate the number of nonzeros in the cholfact factorization of mat, 
along with the number of flops needed to compute it.
Does this through a call to the analyze routine of cholmod.
Note that this is much faster than actually computing the factorization
"""
function ask_cholmod(sdd)
    local (ba, common) = get_suitesparse_handle()

    anal = ba.cholmod_l_analyze(ba.Sparse(sdd), common)
    s_anal = unsafe_load(anal)
    colcount = Base.unsafe_convert(Ptr{Int}, s_anal.ColCount)

    n = Int(s_anal.n)

    nnzL = 0
    flops = 0

    for i = 1:n
        nzl = Int(unsafe_load(colcount, i))
        nnzL += nzl
        flops += nzl^2
    end

    return nnzL, flops
end


"""
  p = cholmod_perm(mat)

Return the permutation that cholmod would apply.
"""
function cholmod_perm(sdd)
    local (ba, common) = get_suitesparse_handle()

    anal = ba.cholmod_l_analyze(ba.Sparse(sdd), common)
    s_anal = unsafe_load(anal)
    perm = Base.unsafe_convert(Ptr{Int}, s_anal.Perm)

    n = Int(s_anal.n)

    p = zeros(Int, n)
    for i = 1:n
        p[i] = unsafe_load(perm, i) + 1
    end

    return p
end
