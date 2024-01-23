#=

approxChol Laplacian solver by Yuan Gao, Rasmus Kyng, and Daniel A. Spielman, 2023.

This algorithm is an implementation of an approximate edge-by-edge elimination
algorithm inspired by the Approximate Gaussian Elimination algorithm of
Kyng and Sachdeva.

A detailed description and experimental evaluation can be found in the paper
"Robust and Practical Solution of Laplacian Equations by Approximate Elimination"
by Yuan Gao, Rasmus Kyng, and Daniel A. Spielman. Paper link: https://arxiv.org/abs/2303.00709.

For usage examples, see https://danspielman.github.io/Laplacians.jl/dev/usingSolvers/

There are several versions of this solver.

As the default solver for Laplacian matrices, we recommend approxchol_lap.
If is not sufficiently reliable for your instances, the more robust approxchol_lap2 can be used.

The solver comes with several different options for orderings, which may
be selected by calling approxchol_lap with different ApproxCholParams(...) options

One version of the solver fixes the order of elimination beforehand,
Another version that adapts the order to eliminate verties of low degree.
These use different data structures.
LLOrdMat is for the fixed order, and LLmatp is for the adaptive order.

These coes produce a structure we call LDLinv that is then used in the solve.
The structure of this code is as follows:

The data structures appear in approxCholTypes.jl
We then have the outline:

* constructors for LLmatp and LLMatOrd
* get_ll_col and compress_ll_col : used inside the elimination
* approxChol : the main routine
* LDLsolver, and its forward and backward solve the apply LDLinv
* approxchol_lap: the main solver, which calls approxchol_lap1 on connected
    components.
    This then calls one of approxchol_lapWdeg, approxchol_lapGiven or approxchol_lapGreedy,
    depending on the parameters.

* approxchol_lapChol - for producing a Cholesky factor instead of an LDLinv.
  might be useful if optimized.
* data structures that are used for the adaptive low-degree version to
  choose the next vertex.

=#

#=
params = ApproxCholParams(order, output, split, merge)

order can be one of
* :deg (by degree, adaptive) -- the default option,
* :wdeg (by original wted degree, nonadaptive),
* :given,
* :random (order by a random permutation)

Furthermore when order option :deg is used,
the two parameters 'split' and 'merge' give addition options.

'split' represents the number of multiedges to split the original edges into.
* If split < 1, it means we are not splitting at all, so we use the initial implementation.
* split = 1 is uses the split-enabled data structures, but does not split the original edges. We recommend not using it.
* split = 2 splits the original edges into 2 copies each. This improves robustness of the solver.

'merge' forces the maximum number of multi-edges
to be bounded above by its value, by merging larger numbers of multi-edges into this number.

It is strongly recommended to use merge == split.

If merge < 1, it means we are not merging at all.
Equivalently, we can set merge to be typemax(Int64),
but that might be slightly less efficient.
Setting merge = 1 will effectively reduce ac 
to its basic form, i.e. split < 1, but setting 
merge to be 1 is less efficient.

=#
mutable struct ApproxCholParams
    order::Symbol
    stag_test::Integer
    split::Integer
    merge::Integer
    keep_last_vertex_for_last::Bool
end

ApproxCholParams() = ApproxCholParams(:deg, 5, 0, 0, false)
ApproxCholParams(sym::Symbol) = ApproxCholParams(sym, 5, 0, 0, false)
ApproxCholParams(sym::Symbol, k) = ApproxCholParams(sym, 5, k, k, false)
ApproxCholParams(sym::Symbol, k, m) = ApproxCholParams(sym, 5, k, m, false)

LDLinv(a::SparseMatrixCSC{Tval,Tind}) where {Tind,Tval} =
  LDLinv(zeros(Tind,a.n-1), zeros(Tind,a.n),Tind[],Tval[],zeros(Tval,a.n))

LDLinv(a::LLMatOrd{Tind,Tval}) where {Tind,Tval} =
  LDLinv(zeros(Tind,a.n-1), zeros(Tind,a.n),Tind[],Tval[],zeros(Tval,a.n))

LDLinv(a::LLmatp{Tind,Tval}) where {Tind,Tval} =
  LDLinv(zeros(Tind,a.n-1), zeros(Tind,a.n),Tind[],Tval[],zeros(Tval,a.n))


function LLmatp(a::SparseMatrixCSC{Tval,Tind}) where {Tind,Tval}
    n = size(a,1)
    m = nnz(a)

    degs = zeros(Tind,n)

    flips = flipIndex(a)

    cols = Array{LLp{Tind,Tval}}(undef, n)
    llelems = Array{LLp{Tind,Tval}}(undef, m)

    @inbounds for i in 1:n
        degs[i] = a.colptr[i+1] - a.colptr[i]

        ind = a.colptr[i]
        j = a.rowval[ind]
        v = a.nzval[ind]
        llpend = LLp{Tind,Tval}(j,v)
        next = llelems[ind] = llpend
        for ind in (a.colptr[i]+one(Tind)):(a.colptr[i+1]-one(Tind))
            j = a.rowval[ind]
            v = a.nzval[ind]
            next = llelems[ind] = LLp{Tind,Tval}(j,v,next)
        end
        cols[i] = next
    end

    @inbounds for i in 1:n
        for ind in a.colptr[i]:(a.colptr[i+1]-one(Tind))
            llelems[ind].reverse = llelems[flips[ind]]
        end
    end

    return LLmatp{Tind,Tval}(n, degs, cols, llelems)
end

function LLmatp(a::SparseMatrixCSC{Tval,Tind}, split::Tind) where {Tind,Tval}
    n = size(a,1)
    m = nnz(a) * split

    degs = zeros(Tind,n)

    flips = flipIndex(a)

    cols = Array{LLp{Tind,Tval}}(undef, n)
    llelems = Array{LLp{Tind,Tval}}(undef, m)

    @inbounds for i in 1:n
        degs[i] = split * (a.colptr[i+1] - a.colptr[i])
        
        # split the edge for split times
        ind = a.colptr[i]
        j = a.rowval[ind]
        v = a.nzval[ind] / split
        llpend = LLp{Tind,Tval}(j,v)
        next = llelems[split * ind] = llpend
        for k in 1:split-1
            next = llelems[split * (ind - 1) + k] = LLp{Tind,Tval}(j,v,next)
        end
        for ind in (a.colptr[i]+one(Tind)):(a.colptr[i+1]-one(Tind))
            
            j = a.rowval[ind]
            v = a.nzval[ind] / split
            # split the edge for split times
            for k in 1:split
                next = llelems[split * (ind - 1) + k] = LLp{Tind,Tval}(j,v,next)
            end
        end
        cols[i] = next
    end

    @inbounds for i in 1:n
        
        for ind in a.colptr[i]:(a.colptr[i+1]-one(Tind))

            # multiedge
            for k in 1:split
                #llelems[split * (ind - 1) + k]
                #llelems[ split * (flips[ind] - 1) + k]
                llelems[split * (ind - 1) + k].reverse = llelems[ split * (flips[ind] - 1) + k]
            end
        end
    end

    return LLmatp{Tind,Tval}(n, degs, cols, llelems)
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

function LLMatOrd(a::SparseMatrixCSC{Tval,Tind}) where {Tind,Tval}
    n = size(a,1)
    m = nnz(a)

    cols = zeros(Tind, n)
    llelems = Array{LLord{Tind,Tval}}(undef, m)

    ptr = one(Tind)

    @inbounds for i in Tind(1):Tind(n-1) # a[n,n] == 0 ?
        next = zero(Tind)

        for ind in (a.colptr[i]):(a.colptr[i+1]-one(Tind))
            j = a.rowval[ind]
            if (i < j)

              v = a.nzval[ind]
              llelems[ptr] = LLord{Tind,Tval}(j, next, v)
              next = ptr
              ptr += one(Tind)

            end
        end
        cols[i] = next
    end

    return LLMatOrd{Tind,Tval}(n, cols, llelems)
end


function LLMatOrd(a::SparseMatrixCSC{Tval,Tind}, split::Tind) where {Tind,Tval}
    n = size(a,1)
    m = nnz(a) * split

    cols = zeros(Tind, n)
    llelems = Array{LLord{Tind,Tval}}(undef, m)

    ptr = one(Tind)

    @inbounds for i in Tind(1):Tind(n-1)
        next = zero(Tind)

        for ind in (a.colptr[i]):(a.colptr[i+1]-one(Tind))
            j = a.rowval[ind]
            if (i < j)
                # split the edge for split times
                v = a.nzval[ind] / split
                for k in 1:split
                    llelems[ptr] = LLord{Tind,Tval}(j, next, v)
                    next = ptr
                    ptr += one(Tind)
                end

            end
        end
        cols[i] = next
    end

    return LLMatOrd{Tind,Tval}(n, cols, llelems)
end


function LLMatOrd(a::SparseMatrixCSC{Tval,Tind}, perm::Array) where {Tind,Tval}
    n = size(a,1)
    m = nnz(a)

    invp = invperm(perm)

    cols = zeros(Tind, n)
    llelems = Array{LLord{Tind,Tval}}(undef, m)

    ptr = one(Tind)

    @inbounds for i0 in Tind(1):Tind(n)
        i = invp[i0]
        next = zero(Tind)

        for ind in (a.colptr[i0]):(a.colptr[i0+1]-one(Tind))
            j = invp[a.rowval[ind]]
            if (i < j)

              v = a.nzval[ind]
              llelems[ptr] = LLord{Tind,Tval}(j, next, v)
              next = ptr
              ptr += one(ptr)

            end
        end
        cols[i] = next
    end

    return LLMatOrd{Tind,Tval}(n, cols, llelems)
end


function LLMatOrd(a::SparseMatrixCSC{Tval,Tind}, perm::Array, split::Tind) where {Tind,Tval}
    n = size(a,1)
    m = nnz(a) * split

    invp = invperm(perm)

    cols = zeros(Tind, n)
    llelems = Array{LLord{Tind,Tval}}(undef, m)

    ptr = one(Tind)

    @inbounds for i0 in Tind(1):Tind(n)
        i = invp[i0]
        next = zero(Tind)

        for ind in (a.colptr[i0]):(a.colptr[i0+1]-one(Tind))
            j = invp[a.rowval[ind]]
            if (i < j)
                # split the edge for split times
                v = a.nzval[ind] / split
                for k in 1:split
                    llelems[ptr] = LLord{Tind,Tval}(j, next, v)
                    next = ptr
                    ptr += one(ptr)
                end

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

function get_ll_col(llmat::LLmatp{Tind,Tval},
  i,
  colspace::Vector{LLp{Tind,Tval}}) where {Tind,Tval}


    ll = llmat.cols[i]
    len = 0
    @inbounds while ll.next != ll

        if ll.val > zero(Tval)
            len = len+1
            if (len > length(colspace))
                push!(colspace,ll)
            else
                colspace[len] = ll
            end
        end

        ll = ll.next
    end

    if ll.val > zero(Tval)
        len = len+1
        if (len > length(colspace))
            push!(colspace,ll)
        else
            colspace[len] = ll
        end
    end

    return len
end

function get_ll_col(llmat::LLMatOrd{Tind,Tval},
  i,
  colspace::Vector{LLcol{Tind,Tval}}) where {Tind,Tval}

    ptr = llmat.cols[i]
    len = 0
    @inbounds while ptr != 0

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



function compressCol!(a::LLmatp{Tind,Tval},
  colspace::Vector{LLp{Tind,Tval}},
  len::Int,
  pq::ApproxCholPQ{Tind}) where {Tind,Tval}

    o = Base.Order.ord(isless, x->x.row, false, Base.Order.Forward)

    sort!(colspace, 1, len, QuickSort, o)

    ptr = 0
    currow::Tind = 0

    c = colspace

    @inbounds for i in 1:len

        if c[i].row != currow
            currow = c[i].row
            ptr = ptr+1
            c[ptr] = c[i]

        else
            c[ptr].val = c[ptr].val + c[i].val
            c[i].reverse.val = zero(Tval)

            approxCholPQDec!(pq, currow)
        end
    end


    o = Base.Order.ord(isless, x->x.val, false, Base.Order.Forward)
    sort!(colspace, 1, ptr, QuickSort, o)

    return ptr
end

function compressCol!(
  colspace::Vector{LLcol{Tind,Tval}},
  len::Int
  ) where {Tind,Tval}

    o = Base.Order.ord(isless, x->x.row, false, Base.Order.Forward)

    sort!(colspace, one(len), len, QuickSort, o)

    c = colspace

    ptr = 0
    currow = c[1].row
    curval = c[1].val
    curptr = c[1].ptr

    @inbounds for i in 2:len

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
    sort!(colspace, one(ptr), ptr, QuickSort, o)

    return ptr
end

function avgCol!(colspace::Vector{LLp{Tind,Tval}},
    len::Int) where {Tind,Tval}

    o_col = Base.Order.ord(isless, x->x.row, false, Base.Order.Forward)
    sort!(colspace, one(len), len, QuickSort, o_col)

    # average the weights of multiedges
    cur_edge = 1
    next_edge = 1

    while next_edge <= len && colspace[cur_edge].row == colspace[next_edge].row
        next_edge += 1
    end
    while cur_edge <= len
        num_mul = next_edge - cur_edge
        mul_sum = 0.0
        for ii in cur_edge : next_edge - 1
            mul_sum += colspace[ii].val
        end
        # change the edge weights
        new_val = mul_sum / num_mul
        for ii in cur_edge : next_edge - 1
            rev_ii = colspace[ii].reverse
            colspace[ii].val = new_val
            rev_ii.val = new_val
        end
        cur_edge = next_edge
        while next_edge <= len && colspace[cur_edge].row == colspace[next_edge].row
            next_edge += 1
        end
    end

    # Lexicographical sorting, first weights, then row indices
    o_col = Base.Order.ord(isless, x->(x.val, x.row), false, Base.Order.Forward)
    sort!(colspace, one(len), len, QuickSort, o_col)
    return len
end

#=
This function averages the weights of multiedges.
If the number of multiedges is larger than merge,
the function also compresses the multiedges 
and keeps only `merge`-multiedges
=#
function compressAvgCol!(colspace::Vector{LLp{Tind,Tval}},
    len::Int,
    pq::ApproxCholPQ{Tind},
    merge::Tind) where {Tind,Tval}

    o_col = Base.Order.ord(isless, x->x.row, false, Base.Order.Forward)
    sort!(colspace, one(len), len, QuickSort, o_col)

    # average the weights of multiedges
    cur_edge = 1
    next_edge = 1

    # ptr points to the end of previous edge.
    # At the end it should be the length of the compressed colspace
    ptr = 0

    while next_edge <= len && colspace[cur_edge].row == colspace[next_edge].row
        next_edge += 1
    end
    while cur_edge <= len
        num_mul = next_edge - cur_edge
        mul_sum = 0.0
        for ii in cur_edge : next_edge - 1
            mul_sum += colspace[ii].val
        end

        if num_mul <= merge
            # change the edge weights
            new_val = mul_sum / num_mul
            for ii in cur_edge : next_edge - 1
                rev_ii = colspace[ii].reverse
                colspace[ii].val = new_val
                rev_ii.val = new_val
                ptr += 1
                colspace[ptr] = colspace[ii]
            end
        else
            # merge and average the multiedges
            new_val = mul_sum / merge
            for ii in cur_edge : cur_edge + merge - 1#next_edge - 1
                rev_ii = colspace[ii].reverse
                colspace[ii].val = new_val
                rev_ii.val = new_val
                ptr += 1
                colspace[ptr] = colspace[ii]
            end
            currow = colspace[cur_edge].row
            for ii in cur_edge + merge : next_edge -1
                colspace[ii].reverse.val = zero(Tval)
                approxCholPQDec!(pq, currow)
            end
        end
        cur_edge = next_edge
        while next_edge <= len && colspace[cur_edge].row == colspace[next_edge].row
            next_edge += 1
        end
    end

    # Lexicographical sorting, first weights, then row indices
    o_col = Base.Order.ord(isless, x->(x.val, x.row), false, Base.Order.Forward)
    sort!(colspace, one(len), ptr, QuickSort, o_col)
    return ptr
end

#=
function approxChol(a::LLMatOrd{Tind,Tval}) where {Tind,Tval}
    n = a.n

    # need to make custom one without col info later.
    ldli = LDLinv(a)
    ldli_row_ptr = one(Tind)

    d = zeros(Tval,n)

    colspace = Array{LLcol{Tind,Tval}}(undef, n)
    cumspace = Array{Tval}(undef, n)
    #vals = Array(Tval,n) # will be able to delete this

    o = Base.Order.ord(isless, identity, false, Base.Order.Forward)


    for i in Tind(1):Tind(n-1)

        ldli.col[i] = i  # will get rid of this with new data type
        ldli.colptr[i] = ldli_row_ptr

        len = get_ll_col(a, i, colspace)

        # RG21 do not compress the column, just sort the entries
        #len = compressCol!(colspace, len)
        o_col = Base.Order.ord(isless, x->x.row, false, Base.Order.Forward)
        sort!(colspace, one(len), len, QuickSort, o_col)


        csum = zero(Tval)
        for ii in 1:len
            #vals[ii] = colspace[ii].val    # if immut, no need for vals
            csum = csum + colspace[ii].val
            # RG21 
            if ii <= length(cumspace)
                cumspace[ii] = csum
            else
                push!(cumspace, csum)
            end
        end
        wdeg = csum

        colScale = one(Tval)

        # RG21 this the edge that is not a multiedge of joffset
        next_edge = 1
        final_d = 0

        for joffset in 1:(len-1)
            # RG21 find the next edge that is not a multiedge of the current one
            isNewNeighbor = false
            if joffset == next_edge
                isNewNeighbor = true
            end
            while next_edge <= len && colspace[joffset].row == colspace[next_edge].row
                next_edge += 1
            end

            # RG21 Avoid adding selfloop at the second last entry
            if next_edge == len + 1
                while joffset <= len - 1
                    ll = colspace[joffset]
                    final_d += ll.val * colScale
                    joffset += 1
                end
                break
            end

            llcol = colspace[joffset]
            w = llcol.val * colScale
            j = llcol.row

            f = w/(wdeg)

            #vals[joffset] = zero(Tval)

            # kind = Laplacians.blockSample(vals,k=1)[1]
            # RG21 sample should start from the next_edge
            r = rand() * (csum - cumspace[next_edge - 1]) + cumspace[next_edge - 1]
            koff = searchsortedfirst(cumspace,r,one(len),len,o)

            k = colspace[koff].row

            #newEdgeVal = w*(one(Tval)-f)
            newEdgeVal = colspace[joffset].val * (csum - cumspace[next_edge - 1]) / csum

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

            #=
            colScale = colScale*(one(Tval)-f)
            #wdeg = wdeg*(1.0-f)^2
            wdeg = wdeg - 2w + w^2/wdeg
            =#

            # RG21 compress all multi edge in ldli 
            if isNewNeighbor#length(ldli.rowval) == 0 || last(ldli.rowval) != j
                push!(ldli.rowval,j)
                push!(ldli.fval, f)
                ldli_row_ptr = ldli_row_ptr + one(Tind)
            else
                ldli.fval[length(ldli.fval)] = last(ldli.fval) + f
            end
            # RG21 only update these after we processed all multiedges
            if joffset == next_edge - 1
                f = last(ldli.fval)
                colScale = colScale*(one(Tval)-f)
                wdeg = wdeg*(one(Tval)-f)^2
            end
            # push!(ops, IJop(i,j,1-f,f))  # another time suck


        end # for joffset


        llcol = colspace[len]
        final_d += llcol.val * colScale
        j = llcol.row

        push!(ldli.rowval,j)
        push!(ldli.fval, one(Tval))
        ldli_row_ptr = ldli_row_ptr + one(Tind)

        d[i] = final_d

    end # for i


    ldli.colptr[n] = ldli_row_ptr

    ldli.d = d

    return ldli
end
=#


function approxChol(a::LLMatOrd{Tind,Tval}) where {Tind,Tval}
    n = a.n

    # need to make custom one without col info later.
    ldli = LDLinv(a)
    ldli_row_ptr = one(Tind)

    d = zeros(Tval,n)

    colspace = Array{LLcol{Tind,Tval}}(undef, n)
    cumspace = Array{Tval}(undef, n)
    #vals = Array(Tval,n) # will be able to delete this

    o = Base.Order.ord(isless, identity, false, Base.Order.Forward)


    for i in Tind(1):Tind(n-1)

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
            koff = searchsortedfirst(cumspace,r,one(len),len,o)

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
            ldli_row_ptr = ldli_row_ptr + one(Tind)

            # push!(ops, IJop(i,j,1-f,f))  # another time suck


        end # for joffset


        llcol = colspace[len]
        w = llcol.val * colScale
        j = llcol.row

        push!(ldli.rowval,j)
        push!(ldli.fval, one(Tval))
        ldli_row_ptr = ldli_row_ptr + one(Tind)

        d[i] = w

    end # for i


    ldli.colptr[n] = ldli_row_ptr

    ldli.d = d

    return ldli
end

# this one is greedy on the degree - also a big win
function approxChol(a::LLmatp{Tind,Tval}) where {Tind,Tval}
    n = a.n

    ldli = LDLinv(a)
    ldli_row_ptr = one(Tind)

    d = zeros(n)

    pq = ApproxCholPQ(a.degs)

    it = 1

    colspace = Array{LLp{Tind,Tval}}(undef, n)
    cumspace = Array{Tval}(undef, n)
    vals = Array{Tval}(undef, n) # will be able to delete this

    o = Base.Order.ord(isless, identity, false, Base.Order.Forward)

    @inbounds while it < n

        i = approxCholPQPop!(pq)

        ldli.col[it] = i # conversion!
        ldli.colptr[it] = ldli_row_ptr

        it = it + 1

        len = get_ll_col(a, i, colspace)

        len = compressCol!(a, colspace, len, pq)  #3hog

        csum = zero(Tval)
        for ii in 1:len
            vals[ii] = colspace[ii].val
            csum = csum + colspace[ii].val
            cumspace[ii] = csum
        end
        wdeg = csum

        colScale = one(Tval)

        for joffset in 1:(len-1)

            ll = colspace[joffset]
            w = vals[joffset] * colScale
            j = ll.row
            revj = ll.reverse

            f = w/(wdeg)

            vals[joffset] = zero(Tval)

            # kind = Laplacians.blockSample(vals,k=1)[1]
            r = rand() * (csum - cumspace[joffset]) + cumspace[joffset]
            koff = searchsortedfirst(cumspace,r,one(len),len,o)

            k = colspace[koff].row

            approxCholPQInc!(pq, k)

            newEdgeVal = f*(one(Tval)-f)*wdeg

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


            colScale = colScale*(one(Tval)-f)
            wdeg = wdeg*(one(Tval)-f)^2

            push!(ldli.rowval,j)
            push!(ldli.fval, f)
            ldli_row_ptr = ldli_row_ptr + one(Tind)

            # push!(ops, IJop(i,j,1-f,f))  # another time suck


        end # for


        ll = colspace[len]
        w = vals[len] * colScale
        j = ll.row
        revj = ll.reverse

        if it < n
            approxCholPQDec!(pq, j)
        end

        revj.val = zero(Tval)

        push!(ldli.rowval,j)
        push!(ldli.fval, one(Tval))
        ldli_row_ptr = ldli_row_ptr + one(Tind)

        d[i] = w

    end

    ldli.colptr[it] = ldli_row_ptr

    ldli.d = d

    return ldli
end

function approxChol(a::LLmatp{Tind,Tval}, split::Int, merge::Int, keep_last_vertex_for_last::Bool = false) where {Tind,Tval}
    n = a.n

    ldli = LDLinv(a)
    ldli_row_ptr = one(Tind)

    d = zeros(n)
    #@show a.degs

    pq = ApproxCholPQ(a.degs, split, keep_last_vertex_for_last)

    it = 1

    colspace = Array{LLp{Tind,Tval}}(undef, n)
    # we use reversed csumspace to avoid subtraction
    csumspace_rev = Array{Tval}(undef, n)
    #cumspace = Array{Tval}(undef, n)
    vals = Array{Tval}(undef, n) # will be able to delete this

    # we use reversed csumspace to avoid subtraction
    #o = Base.Order.ord(isless, identity, false, Base.Order.Forward)
    o_rev = Base.Order.ord(isless, identity, true, Base.Order.Forward)

    @inbounds while it < n

        i = approxCholPQPop!(pq)

        ldli.col[it] = i # conversion!
        ldli.colptr[it] = ldli_row_ptr

        it = it + 1
        #=@show=# len = get_ll_col(a, i, colspace)
        # do not compress the column, just sort the entries
        #len = compressCol!(a, colspace, len, pq)  #3hog
        len = compressAvgCol!(colspace, len, pq, merge)

        csum = zero(Tval)
        for ii in 1:len
            csum = csum + colspace[ii].val
            if ii <= length(csumspace_rev)
                #=@show=# vals[ii] = colspace[ii].val
            else
                push!(vals, colspace[ii].val)
                push!(csumspace_rev, csum)
            end
        end

        csum_rev = zero(Tval)

        for ii in len:-1:1
            csum_rev += colspace[ii].val
            csumspace_rev[ii] = csum_rev
        end


        #=@show=# wdeg = csum

        colScale = one(Tval)

        # next_edge is the edge that is not a multiedge of joffset
        next_edge = 1
        final_d = 0
        for joffset in 1:(len-1)
            # find the next edge that is not a multiedge of the current one
            isNewNeighbor = false
            if joffset == next_edge
                isNewNeighbor = true
            end
            while next_edge <= len && colspace[joffset].row == colspace[next_edge].row
                next_edge += 1
            end

            # Avoid adding selfloop at the second last entry
            if next_edge == len + 1
                # Remove all the multiedges left
                while joffset <= len - 1
                    ll = colspace[joffset]
                    final_d += vals[joffset] * colScale
                    j = ll.row
                    revj = ll.reverse

                    if it < n
                        approxCholPQDec!(pq, j)
                    end

                    revj.val = zero(Tval)
                    joffset += 1
                end
                break
            end


            ll = colspace[joffset]
            
            w = vals[joffset] * colScale
            j = ll.row
            revj = ll.reverse

            #=@show=# f = w/(wdeg)

            vals[joffset] = zero(Tval)

            # kind = Laplacians.blockSample(vals,k=1)[1]
            
            # sample should start from the next_edge
            #=@show=# r = rand() * csumspace_rev[next_edge]
            
            #=@show=# koff = searchsortedlast(csumspace_rev,r,one(len),len,o_rev)

            k = colspace[koff].row

            approxCholPQInc!(pq, k)

            #newEdgeVal = f*(one(Tval)-f)*wdeg
            #newEdgeVal = colspace[joffset].val * (csum - cumspace[next_edge - 1]) / csum
            newEdgeVal = colspace[joffset].val * csumspace_rev[next_edge] / csum

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

            # only update these after we processed all multiedges
            #colScale = colScale*(one(Tval)-f)
            #wdeg = wdeg*(one(Tval)-f)^2

            # compress all multi edge in ldli 
            if isNewNeighbor#length(ldli.rowval) == 0 || last(ldli.rowval) != j
                push!(ldli.rowval,j)
                push!(ldli.fval, f)
                ldli_row_ptr = ldli_row_ptr + one(Tind)
            else
                ldli.fval[length(ldli.fval)] = last(ldli.fval) + f
            end

            # only update these after we processed all multiedges
            if joffset == next_edge - 1
                f = last(ldli.fval)
                #=@show=# colScale = colScale*(one(Tval)-f)
                #=@show=# wdeg = wdeg*(one(Tval)-f)^2
            end

            # push!(ops, IJop(i,j,1-f,f))  # another time suck
        end # for
         
        ll = colspace[len]
        final_d += vals[len] * colScale
        j = ll.row
        revj = ll.reverse

        if it < n
            #@show j
            approxCholPQDec!(pq, j)
        end

        revj.val = zero(Tval)
        
        push!(ldli.rowval,j)
        push!(ldli.fval, one(Tval))
        ldli_row_ptr = ldli_row_ptr + one(Tind)
        
        d[i] = final_d

    end

    ldli.colptr[it] = ldli_row_ptr

    ldli.d = d

    return ldli
end


# This code is deprecated. It should only be used for testing
# Split is used to create ApproxCholPQ
function approxChol(a::LLmatp{Tind,Tval}, split::Int, keep_last_vertex_for_last::Bool = false) where {Tind,Tval}
    n = a.n

    ldli = LDLinv(a)
    ldli_row_ptr = one(Tind)

    d = zeros(n)
    #@show a.degs

    pq = ApproxCholPQ(a.degs, split, keep_last_vertex_for_last)

    it = 1

    colspace = Array{LLp{Tind,Tval}}(undef, n)
    # we use reversed csumspace to avoid subtraction
    csumspace_rev = Array{Tval}(undef, n)
    #cumspace = Array{Tval}(undef, n)
    vals = Array{Tval}(undef, n) # will be able to delete this

    # we use reversed csumspace to avoid subtraction
    #o = Base.Order.ord(isless, identity, false, Base.Order.Forward)
    o_rev = Base.Order.ord(isless, identity, true, Base.Order.Forward)

    @inbounds while it < n

        i = approxCholPQPop!(pq)

        ldli.col[it] = i # conversion!
        ldli.colptr[it] = ldli_row_ptr

        it = it + 1
        #println("In approxChol, iteration $(it - 1), column $(i)")
        #=@show=# len = get_ll_col(a, i, colspace)

        len = avgCol!(colspace, len)

        csum = zero(Tval)
        for ii in 1:len
            csum = csum + colspace[ii].val
            if ii <= length(csumspace_rev)
                #=@show=# vals[ii] = colspace[ii].val
            else
                push!(vals, colspace[ii].val)
                push!(csumspace_rev, csum)
            end
        end

        csum_rev = zero(Tval)

        for ii in len:-1:1
            csum_rev += colspace[ii].val
            csumspace_rev[ii] = csum_rev
        end


        #=@show=# wdeg = csum

        colScale = one(Tval)

        # next_edge is the edge that is not a multiedge of joffset
        next_edge = 1
        final_d = 0
        for joffset in 1:(len-1)
            # find the next edge that is not a multiedge of the current one
            isNewNeighbor = false
            if joffset == next_edge
                isNewNeighbor = true
            end
            while next_edge <= len && colspace[joffset].row == colspace[next_edge].row
                next_edge += 1
            end

            # Avoid adding selfloop at the second last entry
            if next_edge == len + 1
                # Remove all the multiedges left
                while joffset <= len - 1
                    ll = colspace[joffset]
                    final_d += vals[joffset] * colScale
                    j = ll.row
                    revj = ll.reverse

                    if it < n
                        approxCholPQDec!(pq, j)
                    end

                    revj.val = zero(Tval)
                    joffset += 1
                end
                break
            end


            ll = colspace[joffset]
            
            w = vals[joffset] * colScale
            j = ll.row
            revj = ll.reverse

            #=@show=# f = w/(wdeg)

            vals[joffset] = zero(Tval)

            # kind = Laplacians.blockSample(vals,k=1)[1]
            
            # sample should start from the next_edge
            #=@show=# r = rand() * csumspace_rev[next_edge]
            
            #=@show=# koff = searchsortedlast(csumspace_rev,r,one(len),len,o_rev)

            k = colspace[koff].row

            approxCholPQInc!(pq, k)

            #newEdgeVal = f*(one(Tval)-f)*wdeg
            #newEdgeVal = colspace[joffset].val * (csum - cumspace[next_edge - 1]) / csum
            newEdgeVal = colspace[joffset].val * csumspace_rev[next_edge] / csum

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


            # compress all multi edge in ldli 
            if isNewNeighbor#length(ldli.rowval) == 0 || last(ldli.rowval) != j
                push!(ldli.rowval,j)
                push!(ldli.fval, f)
                ldli_row_ptr = ldli_row_ptr + one(Tind)
            else
                ldli.fval[length(ldli.fval)] = last(ldli.fval) + f
            end

            # only update these after we processed all multiedges
            if joffset == next_edge - 1
                f = last(ldli.fval)
                #=@show=# colScale = colScale*(one(Tval)-f)
                #=@show=# wdeg = wdeg*(one(Tval)-f)^2
            end

            # push!(ops, IJop(i,j,1-f,f))  # another time suck
        end # for
        
        ll = colspace[len]
        final_d += vals[len] * colScale
        j = ll.row
        revj = ll.reverse

        if it < n
            #@show j
            approxCholPQDec!(pq, j)
        end

        revj.val = zero(Tval)
        
        push!(ldli.rowval,j)
        push!(ldli.fval, one(Tval))
        ldli_row_ptr = ldli_row_ptr + one(Tind)
        
        d[i] = final_d

    end

    ldli.colptr[it] = ldli_row_ptr

    ldli.d = d

    return ldli
end



#=============================================================

The routines that do the solve.

=============================================================#

function LDLsolver(ldli::LDLinv, b::Vector)
    y = copy(b)

    forward!(ldli, y)

    @inbounds for i in 1:(length(ldli.d))
        if ldli.d[i] != 0
            y[i] /= ldli.d[i]
        end
    end

    backward!(ldli, y)

    mu = mean(y)
    @inbounds for i in eachindex(y)
        y[i] = y[i] - mu
    end

    return y
end


function forward!(ldli::LDLinv{Tind,Tval}, y::Vector) where {Tind,Tval}

    @inbounds for ii in 1:length(ldli.col)
        i = ldli.col[ii]

        j0 = ldli.colptr[ii]
        j1 = ldli.colptr[ii+1]-one(Tind)

        yi = y[i]

        for jj in j0:(j1-1)
            j = ldli.rowval[jj]
            y[j] += ldli.fval[jj] * yi
            yi *= (one(Tval)-ldli.fval[jj])
        end
        j = ldli.rowval[j1]
        y[j] += yi
        y[i] = yi
    end
end

function backward!(ldli::LDLinv{Tind,Tval}, y::Vector) where {Tind,Tval}
    o = one(Tind)
    @inbounds for ii in length(ldli.col):-1:1
        i = ldli.col[ii]

        j0 = ldli.colptr[ii]
        j1 = ldli.colptr[ii+1]-o

        j = ldli.rowval[j1]
        yi = y[i]
        yi = yi + y[j]

        for jj in (j1-o):-o:j0
            j = ldli.rowval[jj]
            yi = (one(Tval)-ldli.fval[jj])*yi + ldli.fval[jj]*y[j]
        end
        y[i] = yi
    end
end

#=
  An attempt at an efficient solver for the case when y is a matrix.
  Have not yet found a meaningful speedup

function LDLsolver(ldli::LDLinv, b::Matrix)
    y = copy(b)

    (d, n) = size(y)
    @assert n == length(ldli.col)+1

    forward!(ldli, y)

    @inbounds for i in 1:(length(ldli.d))
        if ldli.d[i] != 0
          @simd for j in 1:d
            y[j,i] = y[j,i] / ldli.d[i]
          end
        end
    end

    backward!(ldli, y)

    @inbounds for j in 1:size(y,1)
        mu = mean(y[j,:])

        for i in 1:size(y,2)
            y[j,i] = y[j,i] - mu
        end
    end

    return y
end



function forward!{Tind,Tval}(ldli::LDLinv{Tind,Tval}, y::Matrix)

    (d, n) = size(y)
    @assert n == length(ldli.col)+1

    #yi = zeros(y[:,1])

    @inbounds for ii in 1:length(ldli.col)
        i = ldli.col[ii]

        j0 = ldli.colptr[ii]
        j1 = ldli.colptr[ii+1]-one(Tind)

        for jj in j0:(j1-1)
            j = ldli.rowval[jj]
            @simd for k in 1:d
              y[k,j] = y[k,j] + (ldli.fval[jj] * y[k,i])
              y[k,i] = y[k,i] * (one(Tval)-ldli.fval[jj])
          end
        end
        j = ldli.rowval[j1]

        @simd for k in 1:d
          y[k,j] = y[k,j] + y[k,i]
        end
    end
end

function backward!{Tind,Tval}(ldli::LDLinv{Tind,Tval}, y::Matrix)
    o = one(Tind)

    (d, n) = size(y)
    @assert n == length(ldli.col)+1

    yi = zeros(y[:,1])

    @inbounds for ii in length(ldli.col):-1:1
        i = ldli.col[ii]

        j0 = ldli.colptr[ii]
        j1 = ldli.colptr[ii+1]-o

        j = ldli.rowval[j1]
        #copy!(yi, y[:,i])

        @simd for k in 1:d
          y[k,i] = y[k,i] + y[k,j]
        end

        for jj in (j1-o):-o:j0
            j = ldli.rowval[jj]
            @simd for k in 1:d
              y[k,i] = (one(Tval)-ldli.fval[jj])*y[k,i] + ldli.fval[jj].*y[k,j]
            end
        end
        #y[:,i] = yi
    end
end

=#


"""
    solver = approxchol_lap(a); x = solver(b);
    solver = approxchol_lap(a; tol::Real=1e-6, maxits=1000, maxtime=Inf, verbose=false, pcgIts=Int[], params=ApproxCholParams())

A heuristic solver by Yuan Gao, Rasmus Kyng, and Daniel Spielman, see paper https://arxiv.org/abs/2303.00709. The solver is inspired by the solver in https://arxiv.org/abs/1605.02353 by Rasmus Kyng and Sushant Sachdeva. Whereas that paper eliminates vertices one at a time, this eliminates edges one at a time.  It is probably possible to analyze it.
The `ApproxCholParams` let you choose one of three orderings to perform the elimination.

* ApproxCholParams(:given) - in the order given.
    This is the fastest for construction the preconditioner, but the slowest solve.
* ApproxCholParams(:deg) - always eliminate the node of lowest degree.
    This is the slowest build, but the fastest solve.
* ApproxCholParams(:wdeg) - go by a perturbed order of wted degree.

For more info, see http://danspielman.github.io/Laplacians.jl/dev/usingSolvers/index.html
"""
function approxchol_lap(a::SparseMatrixCSC{Tv,Ti};
  tol::Real=1e-6,
  maxits=1000,
  maxtime=Inf,
  verbose=false,
  pcgIts=Int[],
  params=ApproxCholParams()) where {Tv,Ti}

    #@show params.split

    if minimum(a.nzval) < 0
        error("Adjacency matrix can not have negative edge weights")
    end

    return Laplacians.lapWrapComponents(approxchol_lap1, a,
    verbose=verbose,
    tol=tol,
    maxits=maxits,
    maxtime=maxtime,
    pcgIts=pcgIts,
    params=params)


end

"""
    solver = approxchol_lap2(a); x = solver(b);
    solver = approxchol_lap2(a; tol::Real=1e-6, maxits=1000, maxtime=Inf, verbose=false, pcgIts=Int[], params=ApproxCholParams())

`approxchol_lap2` is slower than `approxchol_lap` (by roughly a factor 2), but is more robust.

A heuristic solver by Yuan Gao, Rasmus Kyng, and Daniel Spielman, see paper https://arxiv.org/abs/2303.00709. The solver is inspired by the solver in https://arxiv.org/abs/1605.02353 by Rasmus Kyng and Sushant Sachdeva. Whereas that paper eliminates vertices one at a time, this eliminates edges one at a time.  It is probably possible to analyze it.

`approxchol_lap2` only implements the ordering given by ApproxCholParams(:deg) - always eliminate the node of lowest degree.

For more info, see http://danspielman.github.io/Laplacians.jl/dev/usingSolvers/index.html
"""
function approxchol_lap2(a::SparseMatrixCSC{Tv,Ti};
    tol::Real=1e-6,
    maxits=1000,
    maxtime=Inf,
    verbose=false,
    pcgIts=Int[],
    params=ApproxCholParams()) where {Tv,Ti}
  
      if minimum(a.nzval) < 0
          error("Adjacency matrix can not have negative edge weights")
      end
  
      return Laplacians.lapWrapComponents(approxchol_lap1, a,
      verbose=verbose,
      tol=tol,
      maxits=maxits,
      maxtime=maxtime,
      pcgIts=pcgIts,
      params=ApproxCholParams(:deg, params.stag_test, 2, 2, params.keep_last_vertex_for_last)) # use split = 2 and merge = 2
  
end

function approxchol_lapGreedy(a::SparseMatrixCSC;
  tol::Real=1e-6,
  maxits=1000,
  maxtime=Inf,
  verbose=false,
  pcgIts=Int[],
  params=ApproxCholParams())

  tol_ =tol
  maxits_ =maxits
  maxtime_ =maxtime
  verbose_ =verbose
  pcgIts_ =pcgIts

  t1 = time()

  la = lap(a) # a hit !?


  if params.split >= 1 && params.merge < 1
    llmat = LLmatp(a, params.split)
    ldli = approxChol(llmat, params.split)
  elseif params.split >= 1 && params.merge >= 1
    llmat = LLmatp(a, params.split)
    ldli = approxChol(llmat, params.split, params.merge, params.keep_last_vertex_for_last)
  else
    llmat = LLmatp(a)
    ldli = approxChol(llmat)
  end
  
  
  

  F(b) = LDLsolver(ldli, b)

  if verbose
    println("Using greedy degree ordering. Factorization time: ", time()-t1)
    # println("Condition number of la and ldl is with tol equals 1e-4: $(condNumber(a, ldli; tol=1e-4))")
    println("Factorization size: $(2* length(ldli.fval))")
    println("Ratio of operator edges to original edges: ", 2 * length(ldli.fval) / nnz(a))
    println("ratio of max to min diagonal of laplacian : ", maximum(diag(la))/minimum(diag(la)))
  end

  f(b;tol=tol_,maxits=maxits_, maxtime=maxtime_, verbose=verbose_, pcgIts=pcgIts_) = pcg(la, b .- mean(b), F, tol=tol, maxits=maxits, maxtime=maxtime, pcgIts=pcgIts, verbose=verbose, stag_test = params.stag_test)

end

function approxchol_lapGiven(a::SparseMatrixCSC;
  tol::Real=1e-6,
  maxits=1000,
  maxtime=Inf,
  verbose=false,
  pcgIts=Int[],
  params=ApproxCholParams())

  tol_ =tol
  maxits_ =maxits
  maxtime_ =maxtime
  verbose_ =verbose
  pcgIts_ =pcgIts

  t1 = time()

  la = lap(a)

  llmat = LLMatOrd(a)


  #=
  nv = size(a, 1)
  for i in 1:nv
    print_ll_col(llmat, i)
  end
  println("print all cols complete")
  =#

  ldli = approxChol(llmat)

  

  F(b) = LDLsolver(ldli, b)

  if verbose
    println("Using given ordering. Factorization time: ", time()-t1)
    # println("Condition number of la and ldl tol equals 1e-4: $(condNumber(a, ldli; tol=1e-4))")
    println("Factorization size: $(2* length(ldli.fval))")
    println("Ratio of operator edges to original edges: ", 2 * length(ldli.fval) / nnz(a))
    println("ratio of max to min diagonal of laplacian : ", maximum(diag(la))/minimum(diag(la)))
  end


  f(b;tol=tol_,maxits=maxits_, maxtime=maxtime_, verbose=verbose_, pcgIts=pcgIts_) = pcg(la, b .- mean(b), F, tol=tol, maxits=maxits, maxtime=maxtime, pcgIts=pcgIts, verbose=verbose, stag_test = params.stag_test)

end

function approxchol_lapWdeg(a::SparseMatrixCSC;
  tol::Real=1e-6,
  maxits=1000,
  maxtime=Inf,
  verbose=false,
  pcgIts=Int[],
  params=ApproxCholParams())

  tol_ =tol
  maxits_ =maxits
  maxtime_ =maxtime
  verbose_ =verbose
  pcgIts_ =pcgIts

  t1 = time()

  la = lap(a)

  v = vec(sum(a,dims=1))
  v = v .* (1 .+ rand(length(v)))
  p = sortperm(v)


  llmat = LLMatOrd(a,p)


  #=
  nv = size(a, 1)
  for i in 1:nv
    print_ll_col(llmat, i)
  end
  println("print all cols complete")
  =#

  ldli = approxChol(llmat)

  ip = invperm(p)
  ldlip = LDLinv(p[ldli.col], ldli.colptr, p[ldli.rowval], ldli.fval, ldli.d[ip]);

  F = function(b)
    x = zeros(size(b))
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


  f(b;tol=tol_,maxits=maxits_, maxtime=maxtime_, verbose=verbose_, pcgIts=pcgIts_) = pcg(la, b .- mean(b), F, tol=tol, maxits=maxits, maxtime=maxtime, pcgIts=pcgIts, verbose=verbose, stag_test = params.stag_test)

end

function approxchol_lapRandom(a::SparseMatrixCSC;
    tol::Real=1e-6,
    maxits=1000,
    maxtime=Inf,
    verbose=false,
    pcgIts=Int[],
    params=ApproxCholParams())
  
    tol_ =tol
    maxits_ =maxits
    maxtime_ =maxtime
    verbose_ =verbose
    pcgIts_ =pcgIts
  
    t1 = time()
  
    la = lap(a)
  
    # v = rand(size(a,1))
    # p = sortperm(v)
    p = randperm(size(a,1))
  
    llmat = LLMatOrd(a,p)
  
    #=
    nv = size(a, 1)
    for i in 1:nv
      print_ll_col(llmat, i)
    end
    println("print all cols complete")
    =#
  
    ldli = approxChol(llmat)
  
    ip = invperm(p)
    ldlip = LDLinv(p[ldli.col], ldli.colptr, p[ldli.rowval], ldli.fval, ldli.d[ip]);
  
    F = function(b)
      x = zeros(size(b))
      x = LDLsolver(ldlip, b)
      #x[p] = LDLsolver(ldli, b[p])
      return x
    end
  
    if verbose
        println("Using random ordering. Factorization time: ", time()-t1)
        println("Ratio of operator edges to original edges: ", 2 * length(ldli.fval) / nnz(a))
        println("ratio of max to min diagonal of laplacian : ", maximum(diag(la))/minimum(diag(la)))
    end
  
  
    f(b;tol=tol_,maxits=maxits_, maxtime=maxtime_, verbose=verbose_, pcgIts=pcgIts_) = pcg(la, b .- mean(b), F, tol=tol, maxits=maxits, maxtime=maxtime, pcgIts=pcgIts, verbose=verbose, stag_test = params.stag_test)
  
  end

function approxchol_lap1(a::SparseMatrixCSC{Tv,Ti};
  tol::Real=1e-6,
  maxits=1000,
  maxtime=Inf,
  verbose=false,
  pcgIts=Int[],
  params=ApproxCholParams()) where {Tv,Ti}

    tol_ =tol
    maxits_ =maxits
    maxtime_ =maxtime
    verbose_ =verbose
    pcgIts_ =pcgIts


    if params.order == :deg

      return approxchol_lapGreedy(a,
        verbose=verbose,
        tol=tol,
        maxits=maxits,
        maxtime=maxtime,
        pcgIts=pcgIts,
        params=params)


    elseif params.order == :wdeg

      return approxchol_lapWdeg(a,
        verbose=verbose,
        tol=tol,
        maxits=maxits,
        maxtime=maxtime,
        pcgIts=pcgIts,
        params=params)


    elseif params.order == :random 

        return approxchol_lapRandom(a,
            verbose=verbose,
            tol=tol,
            maxits=maxits,
            maxtime=maxtime,
            pcgIts=pcgIts,
            params=params)
        
    
    else
      return approxchol_lapGiven(a,
        verbose=verbose,
        tol=tol,
        maxits=maxits,
        maxtime=maxtime,
        pcgIts=pcgIts,
        params=params)


    end

end

"""
    solver = approxchol_sddm(sddm); x = solver(b);
    solver = approxchol_sddm(sddm; tol=1e-6, maxits=1000, maxtime=Inf, verbose=false, pcgIts=Int[], params=ApproxCholParams())

Solves sddm systems by wrapping approxchol_lap.
Not yet optimized directly for sddm.

For more info, see http://danspielman.github.io/Laplacians.jl/latest/usingSolvers/index.html
"""
approxchol_sddm = sddmWrapLap(approxchol_lap)
#approxchol_sddm = sddmWrapLapComponents(approxchol_lap1)




#===============================

  Checking the condition number

=================================#

"""
    cn = condNumber(a, ldli; verbose=false)

Given an adjacency matrix a and an ldli computed by approxChol,
this computes the condition number.
"""
function condNumber(a, ldli; verbose=false, tol=1e-2)
  la = lap(a)

  # construct the square operator
  g = function(b)

    y = copy(b)

    #=
    mu = mean(y)
    @inbounds for i in eachindex(y)
        y[i] = y[i] - mu
    end
      =#

    @inbounds for i in 1:(length(ldli.d))
        if ldli.d[i] != 0
            y[i] /= (ldli.d[i])^(1/2)
        else
            y[i] = 0
        end
    end

    backward!(ldli, y)

    y = la * y

    forward!(ldli, y)

    @inbounds for i in 1:(length(ldli.d))
        if ldli.d[i] != 0
            y[i] /= (ldli.d[i])^(1/2)
        else
            y[i] = 0
        end
    end

    #=
    mu = mean(y)
    @inbounds for i in eachindex(y)
        y[i] = y[i] - mu
    end
    =#

    return y
  end

  gOp = SqLinOp(true,1.0,size(a,1),g)
  upper = eigs(gOp;nev=1,which=:LM,tol=tol)[1][1]

  g2(b) = upper*b - g(b)
  g2Op = SqLinOp(true,1.0,size(a,1),g2)
  lower = upper - eigs(g2Op;nev=2,which=:LM,tol=tol)[1][2]

  if verbose
      println("lower: ", lower, ", upper: ", upper);
  end

  return upper/lower

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
        scales[i] = prod(1.0 .- ldli.fval[j0:(j1-1)])
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
    return y .- mean(y)
end


"""
This variation of approxChol creates a cholesky factor to do the elimination.
It has not yet been optimized, and does not yet make the cholesky factor lower triangular
"""
function approxchol_lapChol(a::SparseMatrixCSC{Tv,Ti}; tol::Real=1e-6, maxits=1000, maxtime=Inf, verbose=false, pcgIts=Int[]) where {Tv,Ti}

    tol_ =tol
    maxits_ =maxits
    maxtime_ =maxtime
    verbose_ =verbose
    pcgIts_ =pcgIts

    t1 = time()
    llmat = LLmatp(a)

    ldli = approxChol(llmat)

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


    f(b;tol=tol_,maxits=maxits_, maxtime=maxtime_, verbose=verbose_, pcgIts=pcgIts_) = pcg(la, b .- mean(b), F, tol=tol, maxits=maxits, maxtime=maxtime, pcgIts=pcgIts, verbose=verbose)

    return f
end




#=============================================================

ApproxCholPQ
It only implements pop, increment key, and decrement key.
All nodes with degrees 1 through n appear in their own doubly-linked lists.
Nodes of higher degrees are bundled together.

=============================================================#


function keyMap(x, n)
    return x <= n ? x : n + div(x,n)
end


function keyMap(x, k, upper)
    return x <= k ? x : min(upper, k + div(x, k))
end

function ApproxCholPQ(a::Vector{Tind}, keep_last_vertex_for_last::Bool = false) where Tind

    n = length(a)
    if keep_last_vertex_for_last
        n -= 1
    end

    elems = Array{ApproxCholPQElem{Tind}}(undef, n)
    lists = zeros(Tind, 2*n+1)
    minlist = one(n)

    for i in 1:length(a)
        key = a[i]
        head = lists[key]

        if head > zero(Tind)
            elems[i] = ApproxCholPQElem{Tind}(zero(Tind), head, key)

            elems[head] = ApproxCholPQElem{Tind}(i, elems[head].next, elems[head].key)
        else
            elems[i] = ApproxCholPQElem{Tind}(zero(Tind), zero(Tind), key)

        end

        lists[key] = i
    end

    return ApproxCholPQ(elems, lists, minlist, n, n, one(Tind), keep_last_vertex_for_last ? length(a) : nothing)
end


function ApproxCholPQ(a::Vector{Tind}, split::Int, keep_last_vertex_for_last::Bool = false) where Tind

    n = length(a)
    if keep_last_vertex_for_last
        n -= 1
    end
    elems = Array{ApproxCholPQElem{Tind}}(undef, n)
    lists = zeros(Tind, 2* split * n+1)
    minlist = one(n)

    for i in 1:n
        key = a[i]
        head = lists[key]

        if head > zero(Tind)
            elems[i] = ApproxCholPQElem{Tind}(zero(Tind), head, key)

            elems[head] = ApproxCholPQElem{Tind}(i, elems[head].next, elems[head].key)
        else
            elems[i] = ApproxCholPQElem{Tind}(zero(Tind), zero(Tind), key)

        end

        lists[key] = i
    end

    return ApproxCholPQ(elems, lists, minlist, n, n, split, keep_last_vertex_for_last ? length(a) : nothing)
end

function approxCholPQPop!(pq::ApproxCholPQ{Tind}) where Tind
    if pq.nitems == 0 && !isnothing(pq.last_vertex)
        local v = pq.last_vertex
        pq.last_vertex = nothing
        return v
    end

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
        pq.elems[next] = ApproxCholPQElem(zero(Tind), pq.elems[next].next, pq.elems[next].key)
    end

    pq.nitems -= 1

    return i
end

function approxCholPQMove!(pq::ApproxCholPQ{Tind}, i, newkey, oldlist, newlist) where Tind

    prev = pq.elems[i].prev
    next = pq.elems[i].next

    # remove i from its old list
    if next > zero(Tind)
        pq.elems[next] = ApproxCholPQElem{Tind}(prev, pq.elems[next].next, pq.elems[next].key)
    end
    if prev > zero(Tind)
        pq.elems[prev] = ApproxCholPQElem{Tind}(pq.elems[prev].prev, next, pq.elems[prev].key)

    else
        pq.lists[oldlist] = next
    end

    # insert i into its new list
    head = pq.lists[newlist]
    if head > 0
        pq.elems[head] = ApproxCholPQElem{Tind}(i, pq.elems[head].next, pq.elems[head].key)
    end
    pq.lists[newlist] = i

    pq.elems[i] = ApproxCholPQElem{Tind}(zero(Tind), head, newkey)

    return nothing
end

"""
    Decrement the key of element i
    This could crash if i exceeds the maxkey
"""
function approxCholPQDec!(pq::ApproxCholPQ{Tind}, i) where Tind
    if i == pq.last_vertex
        return nothing
    end

    oldlist = keyMap(pq.elems[i].key, pq.split * pq.n, 2 * pq.split * pq.n + 1)
    newlist = keyMap(pq.elems[i].key - one(Tind), pq.split * pq.n, 2 * pq.split * pq.n + 1)

    if newlist != oldlist

        approxCholPQMove!(pq, i, pq.elems[i].key - one(Tind), oldlist, newlist)

        if newlist < pq.minlist
            pq.minlist = newlist
        end

    else
        pq.elems[i] = ApproxCholPQElem{Tind}(pq.elems[i].prev, pq.elems[i].next, pq.elems[i].key - one(Tind))
    end


    return nothing
end

"""
    Increment the key of element i
    This could crash if i exceeds the maxkey
"""
function approxCholPQInc!(pq::ApproxCholPQ{Tind}, i) where Tind
    if i == pq.last_vertex
        return nothing
    end

    oldlist = keyMap(pq.elems[i].key, pq.split * pq.n, 2 * pq.split * pq.n + 1)
    newlist = keyMap(pq.elems[i].key + one(Tind), pq.split * pq.n, 2 * pq.split * pq.n + 1)

    if newlist != oldlist

        approxCholPQMove!(pq, i, pq.elems[i].key + one(Tind), oldlist, newlist)

    else
        pq.elems[i] = ApproxCholPQElem{Tind}(pq.elems[i].prev, pq.elems[i].next, pq.elems[i].key + one(Tind))
    end

    return nothing
end


