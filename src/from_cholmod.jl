#=
  Routines for finding out what cholmod would do when Cholesky factorizing a matrix

  These are calling cholmod's analyze routinem which is overkill.
  They could be made faster with more direct access to Cholmod.

  By Daniel Spielman.
=#


"""
  nnzL, flops = ask_cholmod(mat)

Estimate the number of nonzeros in the cholfact factorization of mat, 
along with the number of flops needed to compute it.
Does this through a call to the analyze routine of cholmod.
Note that this is much faster than actually computing the factorization
"""
function ask_cholmod(sdd)
    ba = SuiteSparse.CHOLMOD
    #cm = ba.common()
    #cm = ba.defaults(ba.common_struct[Threads.threadid()])

    #anal = ba.analyze(ba.Sparse(lap(sdd)), cm);

    # s_anal = unsafe_load(get(anal.p))
    #s_anal = unsafe_load(pointer(anal))

    anal = ba.cholmod_l_analyze(ba.Sparse(sdd), ba.COMMONS[Threads.threadid()])
    s_anal = unsafe_load(anal)

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
    

"""
  p = cholmod_perm(mat)

Return the permutation that cholmod would apply.
"""
function cholmod_perm(sdd)
    ba = SuiteSparse.CHOLMOD
    #cm = ba.common()
    cm = ba.defaults(ba.common_struct[Threads.threadid()])

    anal = ba.analyze(ba.Sparse(lap(sdd)), cm);

    #s_anal = unsafe_load(get(anal.p))
    s_anal = unsafe_load(pointer(anal))

    n = Int(s_anal.n)

    p = zeros(Int,n)
    for i in 1:n
        p[i] = unsafe_load(s_anal.Perm,i)+1
    end

    return p
end
    
