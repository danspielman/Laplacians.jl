#=

edgeElimOrd Laplacian solver by Daniel A. Spielman, 2017.
This algorithm is an implementation of an approximate edge-by-edge elimination
algorithm inspired by the Approximate Gaussian Elimination algorithm of
Kyng and Sachdeva.

This is for the ordered version, where we know the order before beginning.

=#


LDLinv(n) = LDLinv(zeros(Int,n-1),zeros(Int,n),Array(Int,0),Array(Float64,0),zeros(Float64,n))


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

# note: this code is identical to other, so could probably just strip
# off the types in the function definition
#function get_ll_col(llmat::LLMatOrd, i::Int, colspace::Array{LLord,1})
#
function get_ll_col(llmat, i::Int, colspace)


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


function compressCol!{Tind,Tval}(a::LLMatOrd{Tind,Tval}, colspace::Array{LLord{Tind,Tval},1}, len::Int)

    o = Base.Order.ord(isless, x->x.row, false, Base.Order.Forward)

    sort!(colspace, 1, len, QuickSort, o)

    ptr::Int = 0
    currow::Int = 0

    c = colspace

    for i in 1:len

        # note: if we wanted to implement this with immutables,
        # we could just accumulate the vals and then create object at end

        if c[i].row != currow
            currow = c[i].row
            ptr = ptr+1
            c[ptr] = c[i]

        else
            c[ptr].val = c[ptr].val + c[i].val

        end
    end


    o = Base.Order.ord(isless, x->x.val, false, Base.Order.Forward)
    sort!(colspace, 1, ptr, QuickSort, o)

    return ptr
end


function edgeElim{Tind,Tval}(a::LLMatOrd{Tind,Tval})
    n = a.n

    # need to make custom one without col info later.
    ldli = LDLinv(n)
    ldli_row_ptr = 1

    d = zeros(Tval,n)

    colspace = Array(LLord{Tind,Tval},n)
    cumspace = Array(Tval,n)
    vals = Array(Tval,n) # will be able to delete this

    o = Base.Order.ord(isless, identity, false, Base.Order.Forward)


    # issue: need to allow degree 0.  So, do not want a head,
    # and do not want ptr based.
    for i in 1:(n-1)

        ldli.col[i] = i  # will get rid of this with new data type
        ldli.colptr[i] = ldli_row_ptr

        len = get_ll_col(a, i, colspace)

        len = compressCol!(a,colspace, len)

        csum = 0.0
        for ii in 1:len
            vals[ii] = colspace[ii].val    # if immut, no need for vals
            csum = csum + colspace[ii].val
            cumspace[ii] = csum
        end
        wdeg = csum

        colScale = 1.0

        for joffset in 1:(len-1)

            ll = colspace[joffset]
            w = vals[joffset] * colScale
            j = ll.row

            f = w/(wdeg)

            vals[joffset] = 0.0

            # kind = Laplacians.blockSample(vals,k=1)[1]
            r = rand() * (csum - cumspace[joffset]) + cumspace[joffset]
            koff = searchsortedfirst(cumspace,r,1,len,o)

            k = colspace[koff].row

            newEdgeVal = w*(1.0-f)

            # create edge (j,k) with newEdgeVal
            # do it by reassigning ll
            if j < k # put it in col j
                jhead = a.cols[j]
                ll.next = jhead
                ll.val = newEdgeVal
                ll.row = k
                a.cols[j] = ll
            else # put it in col k
              khead = a.cols[k]
              ll.next = khead
              ll.val = newEdgeVal
              ll.row = j
              a.cols[k] = ll
            end

            colScale = colScale*(1.0-f)
            #wdeg = wdeg*(1.0-f)^2
            wdeg = wdeg - 2w + w^2/wdeg

            push!(ldli.rowval,j)
            push!(ldli.fval, f)
            ldli_row_ptr = ldli_row_ptr + 1

            # push!(ops, IJop(i,j,1-f,f))  # another time suck


        end # for joffset


        ll = colspace[len]
        w = vals[len] * colScale
        j = ll.row

        push!(ldli.rowval,j)
        push!(ldli.fval, 1.0)
        ldli_row_ptr = ldli_row_ptr + 1

        d[i] = w

    end # for i

    ldli.colptr[n] = ldli_row_ptr

    ldli.d = d

    return ldli
end
