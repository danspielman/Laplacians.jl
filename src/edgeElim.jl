#=

edgeElim Laplacian solver by Daniel A. Spielman, 2017.
This algorithm is an implementation of an approximate edge-by-edge elimination
algorithm inspired by the Approximate Gaussian Elimination algorithm of
Kyng and Sachdeva.

There are two versions of this solver:
one that fixes the order of elimination beforehand,
and one that adapts the order to eliminate verties of low degree.
These use different data structures.
LLOrdMat is for the fixed order, and LLmatp is for the adaptive order.

These coes produce a structure we call LDLinv that is then used in the solve.
The structure of this code is as follows:

The data structures appear in EdgeElimTypes.jl
We then have the outline:

* constructors for LLmatp and LLMatOrd
* get_ll_col and compress_ll_col : used inside the elimination
* edgeElim : the main routine
* LDLsolver, and its forward and backward solve the apply LDLinv
* edgeElimLap: the main solver, which calls edgeElimLap1 on connected
    components.
    This then calls one of edgeElimLapWdeg, edgeElimLapGiven or edgeElimLapGreedy,
    depending on the parameters.

* edgeElimLapChol - for producing a Cholesky factor instead of an LDLinv.
  might be useful if optimized.
* data structures that are used for the adaptive low-degree version to
  choose the next vertex.

=#

"""
    params = EdgeElimParams(order, output)
order can be one of
* :deg (by degree, adaptive),
* :wdeg (by original wted degree, nonadaptive),
* :given
"""
type EdgeElimParams
    order::Symbol
end

EdgeElimParams() = EdgeElimParams(:deg)


LDLinv(n) = LDLinv(zeros(Int,n-1),zeros(Int,n),Array(Int,0),Array(Float64,0),zeros(Float64,n))


function LLmatp(a::SparseMatrixCSC)
    n = size(a,1)
    m = nnz(a)

    degs = zeros(Int,n)

    flips = flipIndex(a)

    cols = Array(LLp, n)
    llelems = Array(LLp, m)

    for i in 1:n
        degs[i] = a.colptr[i+1] - a.colptr[i]

        ind = a.colptr[i]
        j = a.rowval[ind]
        v = a.nzval[ind]
        llpend = LLp(j,v)
        next = llelems[ind] = llpend
        for ind in (a.colptr[i]+1):(a.colptr[i+1]-1)
            j = a.rowval[ind]
            v = a.nzval[ind]
            next = llelems[ind] = LLp(j,v,next)
        end
        cols[i] = next
    end

    for i in 1:n
        for ind in a.colptr[i]:(a.colptr[i+1]-1)
            llelems[ind].reverse = llelems[flips[ind]]
        end
    end

    return LLmatp(n, degs, cols, llelems)
end

"""
  Print a column in an LLmatp matrix.
  This is here for diagnostics.
"""
function print_ll_col(llmat::LLmatp, i::Int)
    ll = llmat.cols[i]
    println("col $i, row $(ll.row) : $(ll.val)")

    while ll.next != ll
        ll = ll.next
        println("col $i, row $(ll.row) : $(ll.val)")
    end
end

function LLMatOrd{Tind,Tval}(a::SparseMatrixCSC{Tval,Tind})
    n = size(a,1)
    m = nnz(a)

    cols = zeros(Tind, n)
    llelems = Array(LLord{Tind,Tval}, m)

    ptr = one(Tind)

    for i in 1:(n-1)
        next = zero(Tind)

        for ind in (a.colptr[i]):(a.colptr[i+1]-1)
            j = a.rowval[ind]
            if (i < j)

              v = a.nzval[ind]
              llelems[ptr] = LLord{Tind,Tval}(j, next, v)
              next = ptr
              ptr += 1

            end
        end
        cols[i] = next
    end

    return LLMatOrd{Tind,Tval}(n, cols, llelems)
end

function LLMatOrd{Tind,Tval}(a::SparseMatrixCSC{Tval,Tind}, perm::Array)
    n = size(a,1)
    m = nnz(a)

    invp = invperm(perm)

    cols = zeros(Tind, n)
    llelems = Array(LLord{Tind,Tval}, m)

    ptr = one(Tind)

    for i0 in 1:n
        i = invp[i0]
        next = zero(Tind)

        for ind in (a.colptr[i0]):(a.colptr[i0+1]-1)
            j = invp[a.rowval[ind]]
            if (i < j)

              v = a.nzval[ind]
              llelems[ptr] = LLord{Tind,Tval}(j, next, v)
              next = ptr
              ptr += 1

            end
        end
        cols[i] = next
    end

    return LLMatOrd{Tind,Tval}(n, cols, llelems)
end

"""
  Print a column in an LLMatOrd matrix.
  This is here for diagnostics.
"""
function print_ll_col(llmat::LLMatOrd, i::Int)
    ptr = llmat.cols[i]
    while ptr != 0
      ll = llmat.lles[ptr]
      println("col $i, row $(ll.row) : $(ll.val)")

      ptr = ll.next
    end
end



#=============================================================

The approximate factorization

=============================================================#

function get_ll_col(llmat::LLmatp, i::Int, colspace::Array{LLp,1})


    ll = llmat.cols[i]
    len = 0
    while ll.next != ll

        if ll.val > 0
            len = len+1
            if (len > length(colspace))
                push!(colspace,ll)
            else
                colspace[len] = ll
            end
        end

        ll = ll.next
    end

    if ll.val > 0
        len = len+1
        if (len > length(colspace))
            push!(colspace,ll)
        else
            colspace[len] = ll
        end
    end

    return len
end

function get_ll_col{Tind,Tval}(llmat::LLMatOrd{Tind,Tval}, i::Tind, colspace)

    ptr = llmat.cols[i]
    len = zero(Tind)
    while ptr != 0

        #if ll.val > 0
            len = len+1

            # should not be an lles - is an abuse
            item = LLcol(llmat.lles[ptr].row, ptr, llmat.lles[ptr].val)
            if (len > length(colspace))
                push!(colspace,item)
            else
                colspace[len] = item
            end
        #end

        ptr = llmat.lles[ptr].next
    end

    return len
end



function compressCol!(a::LLmatp, colspace::Array{LLp,1}, len::Int, pq::EdgeElimPQ)

    o = Base.Order.ord(isless, x->x.row, false, Base.Order.Forward)

    sort!(colspace, 1, len, QuickSort, o)

    ptr::Int = 0
    currow::Int = 0

    c = colspace

    for i in 1:len

        if c[i].row != currow
            currow = c[i].row
            ptr = ptr+1
            c[ptr] = c[i]

        else
            c[ptr].val = c[ptr].val + c[i].val
            c[i].reverse.val = 0.0

            edgeElimPQDec!(pq, currow)
        end
    end


    o = Base.Order.ord(isless, x->x.val, false, Base.Order.Forward)
    sort!(colspace, 1, ptr, QuickSort, o)

    return ptr
end

function compressCol!{Tind,Tval}(colspace::Array{LLcol{Tind,Tval},1}, len::Tind)

    o = Base.Order.ord(isless, x->x.row, false, Base.Order.Forward)

    sort!(colspace, 1, len, QuickSort, o)

    c = colspace

    ptr::Tind = 0
    currow = c[1].row
    curval = c[1].val
    curptr = c[1].ptr

    for i in 2:len

        if c[i].row != currow

            ptr = ptr+1
            c[ptr] = LLcol(currow, curptr, curval)  # next is abuse here: reall keep where it came from.

            currow = c[i].row
            curval = c[i].val
            curptr = c[i].ptr

        else

            curval = curval + c[i].val

        end

    end

    # emit the last row

    ptr = ptr+1
    c[ptr] = LLcol(currow, curptr, curval)

    o = Base.Order.ord(isless, x->x.val, false, Base.Order.Forward)
    sort!(colspace, 1, ptr, QuickSort, o)

    return ptr
end


function edgeElim{Tind,Tval}(a::LLMatOrd{Tind,Tval})
    n = a.n

    # need to make custom one without col info later.
    ldli = LDLinv(n)
    ldli_row_ptr = one(Tind)

    d = zeros(Tval,n)

    colspace = Array(LLcol{Tind,Tval},n)
    cumspace = Array(Tval,n)
    #vals = Array(Tval,n) # will be able to delete this

    o = Base.Order.ord(isless, identity, false, Base.Order.Forward)


    for i in 1:(n-1)

        ldli.col[i] = i  # will get rid of this with new data type
        ldli.colptr[i] = ldli_row_ptr

        len = get_ll_col(a, i, colspace)

        len = compressCol!(colspace, len)

        csum = zero(Tval)
        for ii in 1:len
            #vals[ii] = colspace[ii].val    # if immut, no need for vals
            csum = csum + colspace[ii].val
            cumspace[ii] = csum
        end
        wdeg = csum

        colScale = one(Tval)

        for joffset in 1:(len-1)

            llcol = colspace[joffset]
            w = llcol.val * colScale
            j = llcol.row

            f = w/(wdeg)

            #vals[joffset] = zero(Tval)

            # kind = Laplacians.blockSample(vals,k=1)[1]
            r = rand() * (csum - cumspace[joffset]) + cumspace[joffset]
            koff = searchsortedfirst(cumspace,r,1,len,o)

            k = colspace[koff].row

            newEdgeVal = w*(one(Tval)-f)

            # create edge (j,k) with newEdgeVal
            # do it by reassigning ll
            if j < k # put it in col j
                jhead = a.cols[j]
                a.lles[llcol.ptr] = LLord(k, jhead, newEdgeVal)
                #ll.next = jhead
                #ll.val = newEdgeVal
                #ll.row = k
                a.cols[j] = llcol.ptr
            else # put it in col k
              khead = a.cols[k]
              a.lles[llcol.ptr] = LLord(j, khead, newEdgeVal)
              #ll.next = khead
              #ll.val = newEdgeVal
              #ll.row = j
              a.cols[k] = llcol.ptr
            end

            colScale = colScale*(one(Tval)-f)
            #wdeg = wdeg*(1.0-f)^2
            wdeg = wdeg - 2w + w^2/wdeg

            push!(ldli.rowval,j)
            push!(ldli.fval, f)
            ldli_row_ptr = ldli_row_ptr + 1

            # push!(ops, IJop(i,j,1-f,f))  # another time suck


        end # for joffset


        llcol = colspace[len]
        w = llcol.val * colScale
        j = llcol.row

        push!(ldli.rowval,j)
        push!(ldli.fval, one(Tval))
        ldli_row_ptr = ldli_row_ptr + 1

        d[i] = w

    end # for i


    ldli.colptr[n] = ldli_row_ptr

    ldli.d = d

    return ldli
end

# this one is greedy on the degree - also a big win
function edgeElim(a::LLmatp)
    n = a.n

    ldli = LDLinv(n)
    ldli_row_ptr = 1

    d = zeros(n)

    pq = EdgeElimPQ(a.degs)

    it = 1

    colspace = Array(LLp,n)
    cumspace = Array(Float64,n)
    vals = Array(Float64,n) # will be able to delete this

    o = Base.Order.ord(isless, identity, false, Base.Order.Forward)

    while it < n

        i = edgeElimPQPop!(pq)

        ldli.col[it] = i
        ldli.colptr[it] = ldli_row_ptr

        it = it + 1

        len = get_ll_col(a, i, colspace)

        len = compressCol!(a,colspace, len, pq)  #3hog

        csum = 0.0
        for ii in 1:len
            vals[ii] = colspace[ii].val
            csum = csum + colspace[ii].val
            cumspace[ii] = csum
        end
        wdeg = csum

        colScale = 1.0

        for joffset in 1:(len-1)

            ll = colspace[joffset]
            w = vals[joffset] * colScale
            j = ll.row
            revj = ll.reverse

            f = w/(wdeg)

            vals[joffset] = 0.0

            # kind = Laplacians.blockSample(vals,k=1)[1]
            r = rand() * (csum - cumspace[joffset]) + cumspace[joffset]
            koff = searchsortedfirst(cumspace,r,1,len,o)

            k = colspace[koff].row

            edgeElimPQInc!(pq, k)

            newEdgeVal = f*(1-f)*wdeg

            # fix row k in col j
            revj.row = k   # dense time hog: presumably becaus of cache
            revj.val = newEdgeVal
            revj.reverse = ll

            # fix row j in col k
            khead = a.cols[k]
            a.cols[k] = ll
            ll.next = khead
            ll.reverse = revj
            ll.val = newEdgeVal
            ll.row = j


            colScale = colScale*(1-f)
            wdeg = wdeg*(1-f)^2

            push!(ldli.rowval,j)
            push!(ldli.fval, f)
            ldli_row_ptr = ldli_row_ptr + 1

            # push!(ops, IJop(i,j,1-f,f))  # another time suck


        end # for


        ll = colspace[len]
        w = vals[len] * colScale
        j = ll.row
        revj = ll.reverse

        if it < n
            edgeElimPQDec!(pq, j)
        end

        revj.val = 0.0

        push!(ldli.rowval,j)
        push!(ldli.fval, 1.0)
        ldli_row_ptr = ldli_row_ptr + 1

        d[i] = w

    end

    ldli.colptr[it] = ldli_row_ptr

    ldli.d = d

    return ldli
end


#=============================================================

The routines that do the solve.

=============================================================#

function LDLsolver{Tv}(ldli::LDLinv, b::Array{Tv,1})
    y = copy(b)

    forward!(ldli, y)

    for i in 1:(length(ldli.d))
        if ldli.d[i] != 0
            y[i] /= ldli.d[i]
        end
    end

    backward!(ldli, y)

    mu = mean(y)
    for i in eachindex(y)
        y[i] = y[i] - mu
    end

    return y
end


function forward!(ldli::LDLinv, y::Array{Float64,1})
    for ii in 1:length(ldli.col)
        i = ldli.col[ii]

        j0 = ldli.colptr[ii]
        j1 = ldli.colptr[ii+1]-1

        yi = y[i]

        for jj in j0:(j1-1)
            j = ldli.rowval[jj]
            y[j] += ldli.fval[jj] * yi
            yi *= (1-ldli.fval[jj])
        end
        j = ldli.rowval[j1]
        y[j] += yi
        y[i] = yi
    end
end

function backward!(ldli::LDLinv, y::Array{Float64,1})
    for ii in length(ldli.col):-1:1
        i = ldli.col[ii]

        j0 = ldli.colptr[ii]
        j1 = ldli.colptr[ii+1]-1

        j = ldli.rowval[j1]
        yi = y[i]
        yi = yi + y[j]

        for jj in (j1-1):-1:j0
            j = ldli.rowval[jj]
            yi = (1-ldli.fval[jj])*yi + ldli.fval[jj]*y[j]
        end
        y[i] = yi
    end
end





"""
    solver = KMPLapSolver(A; verbose, tol, maxits, maxtime, pcgIts)

A heuristic by Daniel Spielman inspired by the linear system solver in https://arxiv.org/abs/1605.02353 by Rasmus Kyng and Sushant Sachdeva.  Whereas that paper eliminates vertices one at a time, this eliminates edges one at a time.  It is probably possible to analyze it.
The `EdgeElimParams` let you choose one of three orderings to perform the elimination.

* EdgeElimParams(:given) - in the order given.
    This is the fastest for construction the preconditioner, but the slowest solve.
* EdgeElimParams(:deg) - always eliminate the node of lowest degree.
    This is the slowest build, but the fastest solve.
* EdgeElimParams(:wdeg) - go by a perturbed order of wted degree.
    This is the sweet spot in between.
"""
function edgeElimLap{Tv,Ti}(a::SparseMatrixCSC{Tv,Ti};
  tol::Real=1e-6,
  maxits=1000,
  maxtime=Inf,
  verbose=false,
  pcgIts=Int[],
  params=EdgeElimParams())

    return Laplacians.lapWrapComponents(edgeElimLap1, a,
    verbose=verbose,
    tol=tol,
    maxits=maxits,
    maxtime=maxtime,
    pcgIts=pcgIts,
    params=params)


end

function edgeElimLapGreedy(a::SparseMatrixCSC;
  tol::Real=1e-6,
  maxits=1000,
  maxtime=Inf,
  verbose=false,
  pcgIts=Int[],
  params=EdgeElimParams())

  tol_ =tol
  maxits_ =maxits
  maxtime_ =maxtime
  verbose_ =verbose
  pcgIts_ =pcgIts

  t1 = time()

  la = lap(a)

  llmat = LLmatp(a)
  ldli = edgeElim(llmat)
  F(b) = LDLsolver(ldli, b)

  if verbose
    println("Using greedy degree ordering. Factorization time: ", time()-t1)
    println("Ratio of operator edges to original edges: ", 2 * length(ldli.fval) / nnz(a))
  end

  if verbose
      println("ratio of max to min diagonal of laplacian : ", maximum(diag(la))/minimum(diag(la)))
  end


  f(b;tol=tol_,maxits=maxits_, maxtime=maxtime_, verbose=verbose_, pcgIts=pcgIts_) = pcg(la, b-mean(b), F, tol=tol, maxits=maxits, maxtime=maxtime, pcgIts=pcgIts, verbose=verbose)

end

function edgeElimLapGiven(a::SparseMatrixCSC;
  tol::Real=1e-6,
  maxits=1000,
  maxtime=Inf,
  verbose=false,
  pcgIts=Int[],
  params=EdgeElimParams())

  tol_ =tol
  maxits_ =maxits
  maxtime_ =maxtime
  verbose_ =verbose
  pcgIts_ =pcgIts

  t1 = time()

  la = lap(a)

  llmat = LLMatOrd(a)
  ldli = edgeElim(llmat)
  F(b) = LDLsolver(ldli, b)

  if verbose
    println("Using given ordering. Factorization time: ", time()-t1)
    println("Ratio of operator edges to original edges: ", 2 * length(ldli.fval) / nnz(a))
  end

  if verbose
      println("ratio of max to min diagonal of laplacian : ", maximum(diag(la))/minimum(diag(la)))
  end


  f(b;tol=tol_,maxits=maxits_, maxtime=maxtime_, verbose=verbose_, pcgIts=pcgIts_) = pcg(la, b-mean(b), F, tol=tol, maxits=maxits, maxtime=maxtime, pcgIts=pcgIts, verbose=verbose)

end

function edgeElimLapWdeg(a::SparseMatrixCSC;
  tol::Real=1e-6,
  maxits=1000,
  maxtime=Inf,
  verbose=false,
  pcgIts=Int[],
  params=EdgeElimParams())

  tol_ =tol
  maxits_ =maxits
  maxtime_ =maxtime
  verbose_ =verbose
  pcgIts_ =pcgIts

  t1 = time()

  la = lap(a)

  v = vec(sum(a,1))
  v = v .* (1+ rand(length(v)))
  p = sortperm(v)

  llmat = LLMatOrd(a,p)
  ldli = edgeElim(llmat)

  ip = invperm(p)
  ldlip = LDLinv(p[ldli.col], ldli.colptr, p[ldli.rowval], ldli.fval, ldli.d[ip]);

  F = function(b)
    x = zeros(b)
    x = LDLsolver(ldlip, b)
    #x[p] = LDLsolver(ldli, b[p])
    return x
  end

  if verbose
    println("Using wted degree ordering. Factorization time: ", time()-t1)
    println("Ratio of operator edges to original edges: ", 2 * length(ldli.fval) / nnz(a))
  end

  if verbose
      println("ratio of max to min diagonal of laplacian : ", maximum(diag(la))/minimum(diag(la)))
  end


  f(b;tol=tol_,maxits=maxits_, maxtime=maxtime_, verbose=verbose_, pcgIts=pcgIts_) = pcg(la, b-mean(b), F, tol=tol, maxits=maxits, maxtime=maxtime, pcgIts=pcgIts, verbose=verbose)

end



function edgeElimLap1{Tv,Ti}(a::SparseMatrixCSC{Tv,Ti};
  tol::Real=1e-6,
  maxits=1000,
  maxtime=Inf,
  verbose=false,
  pcgIts=Int[],
  params=EdgeElimParams())

    tol_ =tol
    maxits_ =maxits
    maxtime_ =maxtime
    verbose_ =verbose
    pcgIts_ =pcgIts


    if params.order == :deg

      return edgeElimLapGreedy(a,
        verbose=verbose,
        tol=tol,
        maxits=maxits,
        maxtime=maxtime,
        pcgIts=pcgIts,
        params=params)


    elseif params.order == :wdeg

      return edgeElimLapWdeg(a,
        verbose=verbose,
        tol=tol,
        maxits=maxits,
        maxtime=maxtime,
        pcgIts=pcgIts,
        params=params)


    else
      return edgeElimLapGiven(a,
        verbose=verbose,
        tol=tol,
        maxits=maxits,
        maxtime=maxtime,
        pcgIts=pcgIts,
        params=params)


    end

end



#===========================================

  Alternate solver approach

===========================================#


"""
    L = ldli2Chol(ldli)
This produces a matrix L so that L L^T approximate the original Laplacians.
It is not quite a Cholesky factor, because it is off by a perm
(and the all-1s vector orthogonality.
"""
function ldli2Chol(ldli)
    n = length(ldli.colptr)
    m = n + length(ldli.fval)
    li = zeros(Int,m)
    lj = zeros(Int,m)
    lv = zeros(Float64,m)
    lptr = 0

    dhi = zeros(n)
    for i in 1:n
        if ldli.d[i] == 0
            dhi[i] = 1.0
        else
            dhi[i] = sqrt(ldli.d[i])
        end
    end

    scales = ones(n)
    for ii in 1:(n-1)
        i = ldli.col[ii]
        j0 = ldli.colptr[ii]
        j1 = ldli.colptr[ii+1]-1
        scales[i] = prod(1.0-ldli.fval[j0:(j1-1)])
    end

    for ii in 1:(n-1)
        i = ldli.col[ii]
        j0 = ldli.colptr[ii]
        j1 = ldli.colptr[ii+1]-1
        scale = scales[i] / dhi[i]

        scj = 1
        for jj in j0:(j1-1)
            j = ldli.rowval[jj]
            f = ldli.fval[jj]

            lptr += 1
            li[lptr] = i
            lj[lptr] = j
            lv[lptr] = -f*scj/scale


            scj = scj*(1-f)
        end
        j = ldli.rowval[j1]

        lptr += 1
        li[lptr] = i
        lj[lptr] = j
        lv[lptr] = -dhi[i]

        lptr += 1
        li[lptr] = i
        lj[lptr] = i
        lv[lptr] = 1/scale

    end

    for i in 1:n
        if ldli.d[i] == 0
            lptr += 1
            li[lptr] = i
            lj[lptr] = i
            lv[lptr] = 1.0
        end
    end

    return sparse(li,lj,lv,n,n)
    #return li, lj, lv
end

function LDLsolver(L::SparseMatrixCSC, b::Array)
    y = x6 = L \ (L' \ b)
    return y - mean(y)
end


"""
This variation of edgeElim creates a cholesky factor to do the elimination.
It has not yet been optimized, and does not yet make the cholesky factor lower triangular
"""
function edgeElimLapChol{Tv,Ti}(a::SparseMatrixCSC{Tv,Ti}; tol::Real=1e-6, maxits=1000, maxtime=Inf, verbose=false, pcgIts=Int[])

    tol_ =tol
    maxits_ =maxits
    maxtime_ =maxtime
    verbose_ =verbose
    pcgIts_ =pcgIts

    t1 = time()
    llmat = LLmatp(a)

    ldli = edgeElim(llmat)

    chL = ldli2Chol(ldli)

    if verbose
      println("Factorization time: ", time()-t1)
      println("Ratio of operator edges to original edges: ", 2 * length(ldli.fval) / nnz(a))
    end

    F(b) = LDLsolver(chL, b)

    la = lap(a)

    if verbose
        println("ratio of max to min diagonal of laplacian : ", maximum(diag(la))/minimum(diag(la)))
    end


    f(b;tol=tol_,maxits=maxits_, maxtime=maxtime_, verbose=verbose_, pcgIts=pcgIts_) = pcg(la, b-mean(b), F, tol=tol, maxits=maxits, maxtime=maxtime, pcgIts=pcgIts, verbose=verbose)

    return f
end




#=============================================================

EdgeElimPQ
It only implements pop, increment key, and decrement key.
All nodes with degrees 1 through n appear in their own doubly-linked lists.
Nodes of higher degrees are bundled together.

=============================================================#


function keyMap(x::Int, n::Int)
    return x <= n ? x : n + div(x,n)
end

function EdgeElimPQ(a::Array{Int,1})

    n = length(a)
    elems = Array(EdgeElimPQElem,n)
    lists = zeros(Int, 2*n+1)
    minlist = 1

    for i in 1:length(a)
        key = a[i]
        head = lists[key]

        if head > 0
            elems[i] = EdgeElimPQElem(0, head, key)

            elems[head] = EdgeElimPQElem(i, elems[head].next, elems[head].key)
        else
            elems[i] = EdgeElimPQElem(0, 0, key)

        end

        lists[key] = i
    end

    return EdgeElimPQ(elems, lists, minlist, n, n)
end

function edgeElimPQPop!(pq::EdgeElimPQ)
    if pq.nitems == 0
        error("ApproxPQ is empty")
    end
    while pq.lists[pq.minlist] == 0
        pq.minlist = pq.minlist + 1
    end
    i = pq.lists[pq.minlist]
    next = pq.elems[i].next


    pq.lists[pq.minlist] = next
    if next > 0
        pq.elems[next] = EdgeElimPQElem(0, pq.elems[next].next, pq.elems[next].key)
    end

    pq.nitems -= 1

    return i
end

function edgeElimPQMove!(pq::EdgeElimPQ, i::Int, newkey::Int, oldlist::Int, newlist::Int)

    prev = pq.elems[i].prev
    next = pq.elems[i].next

    # remove i from its old list
    if next > 0
        pq.elems[next] = EdgeElimPQElem(prev, pq.elems[next].next, pq.elems[next].key)
    end
    if prev > 0
        pq.elems[prev] = EdgeElimPQElem(pq.elems[prev].prev, next, pq.elems[prev].key)

    else
        pq.lists[oldlist] = next
    end

    # insert i into its new list
    head = pq.lists[newlist]
    if head > 0
        pq.elems[head] = EdgeElimPQElem(i, pq.elems[head].next, pq.elems[head].key)
    end
    pq.lists[newlist] = i

    pq.elems[i] = EdgeElimPQElem(0, head, newkey)

    return Void
end

"""
    Decrement the key of element i
    This could crash if i exceeds the maxkey
"""
function edgeElimPQDec!(pq::EdgeElimPQ, i::Int)

    oldlist = keyMap(pq.elems[i].key, pq.n)
    newlist = keyMap(pq.elems[i].key - 1, pq.n)

    if newlist != oldlist

        edgeElimPQMove!(pq, i, pq.elems[i].key - 1, oldlist, newlist)

        if newlist < pq.minlist
            pq.minlist = newlist
        end

    else
        pq.elems[i] = EdgeElimPQElem(pq.elems[i].prev, pq.elems[i].next, pq.elems[i].key - 1)
    end


    return Void
end

"""
    Increment the key of element i
    This could crash if i exceeds the maxkey
"""
function edgeElimPQInc!(pq::EdgeElimPQ, i::Int)

    oldlist = keyMap(pq.elems[i].key, pq.n)
    newlist = keyMap(pq.elems[i].key + 1, pq.n)

    if newlist != oldlist

        edgeElimPQMove!(pq, i, pq.elems[i].key + 1, oldlist, newlist)

    else
        pq.elems[i] = EdgeElimPQElem(pq.elems[i].prev, pq.elems[i].next, pq.elems[i].key + 1)
    end

    return Void
end
