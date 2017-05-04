#=

edgeElimOrd Laplacian solver by Daniel A. Spielman, 2017.
This algorithm is an implementation of an approximate edge-by-edge elimination
algorithm inspired by the Approximate Gaussian Elimination algorithm of
Kyng and Sachdeva.

This is for the ordered version, where we know the order before beginning.

=#


#LDLinv(n) = LDLinv(zeros(Int,n-1),zeros(Int,n),Array(Int,0),Array(Float64,0),zeros(Float64,n))


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

    invp = zeros(perm)
    invp[perm] = collect(1:n)

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

# note: this code is identical to other, so could probably just strip
# off the types in the function definition
#function get_ll_col(llmat::LLMatOrd, i::Int, colspace::Array{LLcol,1})
#
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

# up to here

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
