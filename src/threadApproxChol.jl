#rked: rasmus editorializing
#TODO: things I need to do
#LATER: things I need to do later


#=

approxChol Laplacian solver by Daniel A. Spielman, 2017.
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

The data structures appear in approxCholTypes.jl
We then have the outline:

* constructors for LLmatp and LLMatOrd
* get_ll_col and compress_ll_col : used inside the elimination
* approxChol : the main routine
* LDLsolver, and its forward and backward solve the apply LDLinv
* approxCholLap: the main solver, which calls approxCholLap1 on connected
components.
This then calls one of approxCholLapWdeg, approxCholLapGiven or approxCholLapGreedy,
depending on the parameters.

* approxCholLapChol - for producing a Cholesky factor instead of an LDLinv.
might be useful if optimized.
* data structures that are used for the adaptive low-degree version to
choose the next vertex.

=#


using DataStructures #TODO #RAT



"""
params = ApproxCholParams(order, output)
order can be one of
* :deg (by degree, adaptive),
* :wdeg (by original wted degree, nonadaptive),
* :given
"""
mutable struct ApproxCholParams
    order::Symbol
    stag_test::Integer
end

ApproxCholParams() = ApproxCholParams(:deg, 5)
ApproxCholParams(sym::Symbol) = ApproxCholParams(sym, 5)

LDLinv{Tind,Tval}(a::SparseMatrixCSC{Tval,Tind}) =
LDLinv(zeros(Tind,a.n-1), zeros(Tind,a.n),Array{Tind}(0),Array{Tval}(0),zeros(Tval,a.n))

LDLinv{Tind,Tval}(a::LLMatOrd{Tind,Tval}) =
LDLinv(zeros(Tind,a.n-1), zeros(Tind,a.n),Array{Tind}(0),Array{Tval}(0),zeros(Tval,a.n))

LDLinv{Tind,Tval}(a::LLmatp{Tind,Tval}) =
LDLinv(zeros(Tind,a.n-1), zeros(Tind,a.n),Array{Tind}(0),Array{Tval}(0),zeros(Tval,a.n))

LDLinv{Tind,Tval}(a::LLmatp{Tind,Tval},m::Tind) =
LDLinv(zeros(Tind,a.n-1), zeros(Tind,a.n),Array{Tind}(m),Array{Tval}(m),zeros(Tval,a.n))

threadLDLinv{Tind,Tval}(a::LLmatp{Tind,Tval},m::Tind) =
threadLDLinv(zeros(Tind,a.n-1), zeros(Tind,a.n-1), zeros(Tind,a.n-1), Array{Tind}(m),Array{Tval}(m),zeros(Tval,a.n))

function LLmatp{Tind,Tval}(a::SparseMatrixCSC{Tval,Tind})
    n = size(a,1)
    m = nnz(a)

    degs = zeros(Tind,n)

    flips = flipIndex(a)

    cols = Array{LLp{Tind,Tval}}(n)
    llelems = Array{LLp{Tind,Tval}}(m)

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

function threadPrintLLcol(llmat::LLmatp, i::Int)
    #ccall(:jl_,Void,(Any,), "pid=$(Threads.threadid()), threadPrintLLcol")
    ll = llmat.cols[i]
    #ccall(:jl_,Void,(Any,), "pid=$(Threads.threadid()), col $i, row $(ll.row) : $(ll.val)")
    while ll.next != ll
        ll = ll.next
        #ccall(:jl_,Void,(Any,), "pid=$(Threads.threadid()), col $i, row $(ll.row) : $(ll.val)")
    end
end

function LLMatOrd{Tind,Tval}(a::SparseMatrixCSC{Tval,Tind})
    n = size(a,1)
    m = nnz(a)

    cols = zeros(Tind, n)
    llelems = Array{LLord{Tind,Tval}}(m)

    ptr = one(Tind)

    @inbounds for i in Tind(1):Tind(n-1)
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

function LLMatOrd{Tind,Tval}(a::SparseMatrixCSC{Tval,Tind}, perm::Array)
    n = size(a,1)
    m = nnz(a)

    invp = invperm(perm)

    cols = zeros(Tind, n)
    llelems = Array{LLord{Tind,Tval}}(m)

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

function claimForElim{Tind}(threadDebugState::ThreadDebugState, pState::Array{Atomic{Tind}},i::Tind,
    claimStack::Array{Tind}, vState::Array{Atomic{Tind}},vLock::Array{Base.Threads.AbstractLock})
tryWhileCount = 0
claimForElimResult = zero(Tind)    
while claimForElimResult == zero(Tind)
    claimForElimResult = claimVertex(threadDebugState,pState,claimStack,vState,vLock,i)
    if claimForElimResult == 0 #retry
        lock(vLock[i])
        unlock(vLock[i])
        #TODO is this lock/unlock sensible?
    end
    if tryWhileCount > 10000
        #ccall(:jl_,Void,(Any,), "pid=$(Threads.threadid()), at MARK cFE 1")
        @assert false "pid=$(Threads.threadid()): too many tries in claimForElim"
    end
    tryWhileCount += 1
end
if claimForElimResult == -one(Tind)
    return false
end
if claimForElimResult == one(Tind)
    return true
end
assert(false) #should not reach here #RAT
return false
end


#rked claim the nbrs
function claimAllNbrs{Tind,Tval}(threadDebugState::ThreadDebugState, pState::Array{Atomic{Tind}},llmat::LLmatp{Tind,Tval},i::Tind,claimStack::Array{Tind}, vState::Array{Atomic{Tind}},vLock::Array{Base.Threads.AbstractLock})

    # at this point, we should be guaranteed that llmat.cols[i] belongs to us,
    # i.e. no one else will modify it
    # so the set of nbrs is well-defined
    tryWhileCount = 0
    ll = llmat.cols[i]
    
    doWhileTest = true
    claimNbrResult =  zero(Tind)
    #TODO @inbounds?
    while doWhileTest
        if tryWhileCount > 10000
            ##ccall(:jl_,Void,(Any,), "pid=$(Threads.threadid()), at MARK cAN 1")
            printDebug(threadDebugState)
            @assert false "pid=$(Threads.threadid()): too many tries in claimAllNbrs"
        end
        tryWhileCount += 1
        if ll.val > zero(Tval)
            v = ll.row
            claimNbrResult = claimVertex(threadDebugState,pState,claimStack,vState,vLock,v)
            if claimNbrResult == zero(Tind) #retry
                #RJK
                ##ccall(:jl_,Void,(Any,), "pid=$(Threads.threadid()), waiting on vLock[$(v)]")
                #lock(vLock[v])
                #unlock(vLock[v]) #WAIT THIS IS PROBABLY NOT OK?
                ##ccall(:jl_,Void,(Any,), "pid=$(Threads.threadid()), released vLock[$(v)]")
                continue #retry
            end
            if claimNbrResult == -one(Tind) #failure
                return false
            end
            @assert claimNbrResult == one(Tind) "claimNbr success" #success
        end
        doWhileTest = (ll.next != ll)
        ll = ll.next
    end

    return true
end

#rked extract copy of col from overall memspace
function get_ll_col{Tind,Tval}(llmat::LLmatp{Tind,Tval},
 i,
 colspace::Vector{LLp{Tind,Tval}})
ll = llmat.cols[i]
len = 0
while ll.next != ll
    #@inbounds while ll.next != ll #RAT back to @inbounds

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

function get_ll_col{Tind,Tval}(llmat::LLMatOrd{Tind,Tval},
 i,
 colspace::Vector{LLcol{Tind,Tval}})

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

function TESTPQcompressCol!{Tind,Tval}(a::LLmatp{Tind,Tval},
   colspace::Vector{LLp{Tind,Tval}}, #rked #TODO Q is this a pointer to inside a larger structure, or is it a disposable chunk?
   len::Int)
# len::Int,
# pq::ApproxCholPQ{Tind},
# pqLock::Base.Threads.AbstractLock)

#ccall(:jl_,Void,(Any,), "pid=$(Threads.threadid()), at TESTPQcompressCol")

o = Base.Order.ord(isless, x->x.row, false, Base.Order.Forward)

sort!(colspace, 1, len, QuickSort, o)

ptr = 0
currow::Tind = 0

c = colspace

for i in 1:len 
    #@inbounds for i in 1:len #RAT back to inbounds

    if c[i].row != currow
        currow = c[i].row
        ptr = ptr+1
        c[ptr] = c[i] #rked ah, we start treating the first part of array as a new array

    else
        c[ptr].val = c[ptr].val + c[i].val
        c[i].reverse.val = zero(Tval) #TODO careful! this relies on us owning the list where reverse lives
        #rked is this how the reverse edge is killed
        #rked #TODO how is the space later recovered?
        #rked #TODO why is c[i].val not zeroed out?

        # #TODO MUTEX NOT OK HERE?!
        # lock(pqLock)
        # approxCholPQDec!(pq, currow)
        # atomic_fence()
        # unlock(pqLock)
    end
end


o = Base.Order.ord(isless, x->x.val, false, Base.Order.Forward) 
sort!(colspace, 1, ptr, QuickSort, o) #rked we then sort only the part that's been allocated for the "new array"?

return ptr #rked looks like it's returning the vector with length noted, and then rest is considered useless
end


#rked Q why 'a' input?
#rked A unclear: doesn't seem to be used but
#rked #LATER confirm
function compressCol!{Tind,Tval}(a::LLmatp{Tind,Tval},
   colspace::Vector{LLp{Tind,Tval}}, #rked #TODO Q is this a pointer to inside a larger structure, or is it a disposable chunk?
   len::Int,
   pq::ApproxCholPQ{Tind},
   pqLock::Base.Threads.AbstractLock)

o = Base.Order.ord(isless, x->x.row, false, Base.Order.Forward)

sort!(colspace, 1, len, QuickSort, o)

ptr = 0
currow::Tind = 0

c = colspace

for i in 1:len 
    #@inbounds for i in 1:len #RAT back to inbounds

    if c[i].row != currow
        currow = c[i].row
        ptr = ptr+1
        c[ptr] = c[i] #rked ah, we start treating the first part of array as a new array

    else
        c[ptr].val = c[ptr].val + c[i].val
        c[i].reverse.val = zero(Tval) #TODO careful! this relies on us owning the list where reverse lives
        #rked is this how the reverse edge is killed
        #rked #TODO how is the space later recovered?
        #rked #TODO why is c[i].val not zeroed out?

        #TODO MUTEX NOT OK HERE?!
        lock(pqLock)
        approxCholPQDec!(pq, currow)
        atomic_fence()
        unlock(pqLock)
    end
end


o = Base.Order.ord(isless, x->x.val, false, Base.Order.Forward) 
sort!(colspace, 1, ptr, QuickSort, o) #rked we then sort only the part that's been allocated for the "new array"?

return ptr #rked looks like it's returning the vector with length noted, and then rest is considered useless
end

function compressCol!{Tind,Tval}(a::LLmatp{Tind,Tval},
   colspace::Vector{LLp{Tind,Tval}}, #rked #TODO Q is this a pointer to inside a larger structure, or is it a disposable chunk?
   len::Int,
   pq::ApproxCholPQ{Tind})

o = Base.Order.ord(isless, x->x.row, false, Base.Order.Forward)

sort!(colspace, 1, len, QuickSort, o)

ptr = 0
currow::Tind = 0

c = colspace

for i in 1:len 
    #@inbounds for i in 1:len #RAT back to inbounds

    if c[i].row != currow
        currow = c[i].row
        ptr = ptr+1
        c[ptr] = c[i] #rked ah, we start treating the first part of array as a new array

    else
        c[ptr].val = c[ptr].val + c[i].val
        c[i].reverse.val = zero(Tval) #TODO careful! this relies on us owning the list where reverse lives
        #rked is this how the reverse edge is killed
        #rked #TODO how is the space later recovered?
        #rked #TODO why is c[i].val not zeroed out?

        approxCholPQDec!(pq, currow)
    end
end


o = Base.Order.ord(isless, x->x.val, false, Base.Order.Forward) 
sort!(colspace, 1, ptr, QuickSort, o) #rked we then sort only the part that's been allocated for the "new array"?

return ptr #rked looks like it's returning the vector with length noted, and then rest is considered useless
end

function compressCol!{Tind,Tval}(
   colspace::Vector{LLcol{Tind,Tval}},
   len::Int
   )

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


function approxChol{Tind,Tval}(a::LLMatOrd{Tind,Tval})
    n = a.n

    # need to make custom one without col info later.
    ldli = LDLinv(a)
    ldli_row_ptr = one(Tind)

    d = zeros(Tval,n)

    colspace = Array{LLcol{Tind,Tval}}(n)
    cumspace = Array{Tval}(n)
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
function approxChol{Tind,Tval}(a::LLmatp{Tind,Tval})
    n = a.n

    ldli = LDLinv(a)
    ldli_row_ptr = one(Tind)

    d = zeros(n)

    pq = ApproxCholPQ(a.degs)

    it = 1

    colspace = Array{LLp{Tind,Tval}}(n)
    cumspace = Array{Tval}(n)
    vals = Array{Tval}(n) # will be able to delete this

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

function testThreadApproxCholLapGreedy(a::SparseMatrixCSC;
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

llmat = LLmatp(a)
ldli = testThreadApproxChol(llmat,nnz(a))
F(b) = threadLDLsolver(ldli, b)

if verbose
    println("Using greedy degree ordering. Factorization time: ", time()-t1)
    println("Ratio of operator edges to original edges: ", 2 * length(ldli.fval) / nnz(a))
end

if verbose
    println("ratio of max to min diagonal of laplacian : ", maximum(diag(la))/minimum(diag(la)))
end


f(b;tol=tol_,maxits=maxits_, maxtime=maxtime_, verbose=verbose_, pcgIts=pcgIts_) = pcg(la, b-mean(b), F, tol=tol, maxits=maxits, maxtime=maxtime, pcgIts=pcgIts, verbose=verbose, stag_test = params.stag_test)

return (f,ldli)

end

function threadApproxCholLapGreedy(a::SparseMatrixCSC;
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

llmat = LLmatp(a)
ldli = threadApproxChol(llmat,nnz(a))
F(b) = LDLsolver(ldli, b)

if verbose
    println("Using greedy degree ordering. Factorization time: ", time()-t1)
    println("Ratio of operator edges to original edges: ", 2 * length(ldli.fval) / nnz(a))
end

if verbose
    println("ratio of max to min diagonal of laplacian : ", maximum(diag(la))/minimum(diag(la)))
end


f(b;tol=tol_,maxits=maxits_, maxtime=maxtime_, verbose=verbose_, pcgIts=pcgIts_) = pcg(la, b-mean(b), F, tol=tol, maxits=maxits, maxtime=maxtime, pcgIts=pcgIts, verbose=verbose, stag_test = params.stag_test)

end

#TODO are these of acceptable types?! probably not
const PS_ASSISTING = convert(Int64,0)
const PS_CLAIMING = convert(Int64,1)
const PS_ELIMINATING = convert(Int64,2)

#TODO careful: we distinguish between "pHasVert[threadid()] == i" and i being claimed...
function claimVertex{Tind}(threadDebugState::ThreadDebugState, pState::Array{Atomic{Tind}}, claimStack::Array{Tind},vState::Array{Atomic{Tind}},
 vLock::Array{Base.Threads.AbstractLock}, v::Tind)
#TODO should we check our own thread state occasionally? e.g. like selfState below
#TODO or should we remove it for speed?
pushDebugMsg!(threadDebugState, "pid=$(Threads.threadid()) begins claim $(v)")
selfState = pState[threadid()][]
assert( selfState != PS_ELIMINATING )
if selfState == PS_ASSISTING
    pushDebugMsg!(threadDebugState, "pid=$(Threads.threadid()) assisting, so failed claim")
    return -one(Tind) #failure
end
vsOld = atomic_cas!(vState[v], zero(Tind), threadid()) #0 is for unclaimed
assert(vsOld != -1) #TODO we don't try to claim eliminated vertices?  
if vsOld == zero(Tind)
    assert(vState[v][] == threadid()) # thread must now have claimed the vertex #RAT
    push!(claimStack,v)
    lock(vLock[v]) #TODO note I believe this might be important for memory consistency?
    #not sure how that compares to atomics --- does a homemade spin lock also give memory consistency?
    pushDebugMsg!(threadDebugState, "pid=$(Threads.threadid()) claimed $(v)")
    return one(Tind) #success
end
if vsOld == threadid() # thread had already previously claimed the vertex #RAT
    assert(vState[v][] == threadid()) #RAT
    pushDebugMsg!(threadDebugState, "pid=$(Threads.threadid()) ALREADY claimed $(v) - could happen if there are duplicate edges, as this is done before dedup")
    #@assert false
    return one(Tind) #success
end
if vsOld < threadid()
    # a thread with higher priority (lower id) has claimed the vertex and we need to abandon
    pushDebugMsg!(threadDebugState, "pid=$(Threads.threadid()) lost priority, so failed claim")
    return -one(Tind) #failure
end
if vsOld > threadid()
    pushDebugMsg!(threadDebugState, "pid=$(Threads.threadid()) is trying to get priority over thread $(vsOld)")
    # a thread with higher priority (lower id) has claimed the vertex
    # we should try to claim it, and fail only if the other has reached the eliminating state
    competingPsOld = atomic_cas!(pState[vsOld],PS_CLAIMING,PS_ASSISTING) #0 is for unclaimed
    if competingPsOld == PS_ELIMINATING
        pushDebugMsg!(threadDebugState, "pid=$(Threads.threadid()) failed to get priority over thread $(vsOld)")
        # the other thread reach the stage of eliminating
        # so we should abandon... #right?
        #TODO or could we wait on the vertex to become accessible?...
        #TODO we know it is not being eliminated, because then it would have been adj to our vertex, which we own...
        #TODO EXCEPT: what if we are still trying to get first vertex?
        pushDebugMsg!(threadDebugState, "pid=$(Threads.threadid()) competitor reached ELIM state, so failed claim")
        return -one(Tind) #failure
    end
    if competingPsOld == PS_CLAIMING || competingPsOld == PS_ASSISTING
        pushDebugMsg!(threadDebugState, "pid=$(Threads.threadid()) got priority over thread $(vsOld)")
        #we should retry
        #TODO but I think there's an issue: what if the vertex was since eliminated?
        #CLAIM but it can't have been?
        #If it is our attempt to claim the vertex we eliminate,
        #then no one else will try to eliminate it, because elim ownership goes via the PQ mutex
        #If it is a DIFFERENT vertex, v != i, then we must already own the vertex we are trying to eliminate (i),
        # and v is a nbr of i, so v cannot be elim'd while we own i
        pushDebugMsg!(threadDebugState, "pid=$(Threads.threadid()) waiting for competitor to assist, so retry claim")
        return zero(Tind) #retry
    end
end

assert(false) #we should never get to this line
return zero(Tind)
end


function releaseStack{Tind}(threadDebugState::ThreadDebugState, vState::Array{Atomic{Tind}},vLock::Array{Base.Threads.AbstractLock},claimStack::Array{Tind})
    for s in eachindex(claimStack)
        v = pop!(claimStack)
        pushDebugMsg!(threadDebugState, "pid=$(Threads.threadid()) releases $(v)")
        vs = vState[v]
        assert(vs[] == threadid())
        vs[] = 0 #UNCLAIMED
        unlock(vLock[v])
    end
end

function myAbandon{Tind}(threadDebugState::ThreadDebugState, pState::Atomic{Tind}, vState::Array{Atomic{Tind}}, vLock::Array{Base.Threads.AbstractLock},
   claimStack::Array{Tind})
pState[] = PS_ASSISTING #TODO do we want that here?
pushDebugMsg!(threadDebugState, "pid=$(Threads.threadid()), abandoning POST assist assert")
releaseStack(threadDebugState, vState,vLock,claimStack)
#NB: pHasVert[threadid()] is till set to i
end


#TODO WARNING: currently I'm assuming Tind allows for negative numbers?
# this one is greedy on the degree - also a big win
function threadApproxChol{Tind,Tval}(a::LLmatp{Tind,Tval},m::Tind)

end

function testThreadApproxChol{Tind,Tval}(a::LLmatp{Tind,Tval},m::Tind)
    pCount = nthreads() #process count

    # state for processors
    pState = Array{Atomic{Tind}}(pCount)
    for p = 1:pCount
        pState[p] = Atomic{Tind}(PS_ASSISTING)
    end
    pHasVert = zeros(Tind,pCount) # 0 for when p as no vertex, v > 0 for vertex id otherwise

    n = a.n
    # state for vertices
    vLock = Array{Base.Threads.AbstractLock}(n)
    vState = Array{Atomic{Int64}}(n) #TODO locks mean we don't need atomic anymore?
    for v = 1:n
        vState[v] = Atomic{Int64}(0) #UNCLAIMED
        vLock[v] = SpinLock()
    end

    #lock for priority queue access (variable pq)
    assert(pCount < 40) # we need to switch from SpinLock if using >> 30 processors?
    pqLock = SpinLock() #THREADGLOBAL     

    ldli = threadLDLinv(a,2*m) #THREADGLOBAL #TODO I don't know if this will work with m or m/2 or 10m !!!
    ldli_row_ptr = Atomic{Tind}(one(Tind)) #TODO #THREADGLOBAL?

    d = zeros(n) #THREADGLOBAL

    pq = ThreadApproxCholPQ(a.degs)  #THREADGLOBAL

    pDebugMsgQueues = Array{CircularBuffer{String}}(pCount)
    for p = 1:pCount
        pDebugMsgQueues[p] = CircularBuffer{String}(100)
    end
    threadDebugState = ThreadDebugState(pDebugMsgQueues,pState,vLock,vState,pq,pHasVert,Atomic{Int64}(1))
    pushDebugMsg!(threadDebugState,"test msg")

    atomicIt = Atomic{Int64}(1) #THREADGLOBAL
    atomicCountThreadsContinuing = Atomic{Int64}(0) #THREADGLOBAL

    #TODO maybe we could save some indirection
    # by allocating these things inside a threaded loop
    # and then running while loops inside that loop for actual elim?
    colspace = Array{Array{LLp{Tind,Tval}}}(pCount) #TODO #THREADLOCAL
    for pid = 1:pCount #TODO could @threads this
        colspace[pid] = Array{LLp{Tind,Tval}}(n)
    end    
    cumspace = Array{Array{Tval}}(pCount)  #TODO #THREADLOCAL
    for pid = 1:pCount #TODO could @threads this
        cumspace[pid] = Array{Tval}(n)
    end 
    vals = Array{Array{Tval}}(pCount)  #TODO #THREADLOCAL
    for pid = 1:pCount #TODO could @threads this
        vals[pid] = Array{Tval}(n)
    end

    claimStack = Array{Array{Tind}}(pCount)  #TODO #THREADLOCAL
    for pid = 1:pCount #TODO could @threads this
        claimStack[pid] = Array{Tind}(0)
    end

    rngSeeds = rand(UInt64,pCount)  #TODO #THREADLOCAL
    
    o = Base.Order.ord(isless, identity, false, Base.Order.Forward) #TODO #THREADGLOBAL? UNSURE, but hopefully?

    #RAT
    # This should kill attempts by proc 1 to claim vertex 3
    # vState[3][] = 2
    # #pState[2][] = PS_ELIMINATING
    # #pState[2][] = PS_CLAIMING
    # pState[2][] = PS_ASSISTING


    #@inbounds
    #maxCollisionCount = 5*n #TODO reduce by a lot
    #for procId = 1:nthreads()
    @threads for procId = 1:nthreads()
    #for procId = 1:1

    pushDebugMsg!(threadDebugState, "pid=$(Threads.threadid()) handles for loop iter procId=$(procId)")
    rng = MersenneTwister(rngSeeds[threadid()])
    it = 0 #TODO #THREADLOCAL
    #TODO temp doing counter this way
    tryWhileThreshold = 2000000
    tryWhileCount = 0 #tryWhileCount
    while tryWhileCount < tryWhileThreshold #200 also didn't hang, at least sometimes
        tryWhileCount += 1
        @assert tryWhileCount < tryWhileThreshold
        #@threads for loopIt = 1:(n-1+maxCollisionCount)
        #for loopIt = 1:(n-1+maxCollisionCount)
        # this logic for it = atomic_add... probably doesn't actually work when multi-threaded!
        # instead, we need to count number of eliminations done
        # if threadid() != 2  #RAT
        #     atomic_add!(atomicCountThreadsContinuing,1)
        #     continue
        # end

        pushDebugMsg!(threadDebugState, "pid=$(Threads.threadid()), starting for loop")

        #pushDebugMsg!(threadDebugState, "pid=$(Threads.threadid()), at MARK 1\n")
        pState[threadid()][] = PS_CLAIMING #THREADGLOBAL 
        i = pHasVert[threadid()]
        pushDebugMsg!(threadDebugState, "pid=$(Threads.threadid()), read hasVert as i=$(i)")
        if i == 0
            it = atomic_add!(atomicIt,1)
            if it > n-1
                atomic_add!(atomicCountThreadsContinuing,1)
                break #HIGHTODO #TODO #RJK WORKING HERE  THINK THIS WILL EVENTUALLY EXIT ALL THREADS
            end

            #MUTEX enter mutual exclusion region
            lock(pqLock)
            pushDebugMsg!(threadDebugState,"pid=$(Threads.threadid()), PRE-pop-print") #RAT
            #printPQ(pq)
            #TODO need to allow for skipping if empty here
            i = threadApproxCholPQPop!(pq) #THREADLOCAL i #localscope #TESTPQ
            pushDebugMsg!(threadDebugState,"pid=$(Threads.threadid()), popped i=$(i)") #RAT
            pushDebugMsg!(threadDebugState,"pid=$(Threads.threadid()), POST-pop-print") #RAT
            #printPQ(pq)
            unlock(pqLock)

            #EXIT mutual exclusion region
            pHasVert[threadid()] = i
        end

        #TODO
        # Here I've used very infrequent checks for state switch
        # would probably be better to at least sometimes check if another proc

        pushDebugMsg!(threadDebugState, "pid=$(Threads.threadid()), at PRE-claimForElim, i=$(i), pHasVert[threadid()]=$(pHasVert[threadid()])")

        # #pushDebugMsg!(threadDebugState, "pid=$(Threads.threadid()), at MARK 2")

        claimForElimSucceeded = claimForElim(threadDebugState, pState,i,claimStack[threadid()],vState,vLock)
        #@assert claimForElimSucceeded "claimForElimSucceeded"  #success #RAT
        if !claimForElimSucceeded #failure
            myAbandon(threadDebugState, pState[threadid()],vState,vLock,claimStack[threadid()])
            continue  #leave the for loop #TODO: careful if switching to top-level for + inside while
        end
        pushDebugMsg!(threadDebugState, "pid=$(Threads.threadid()), at POST-claimForElim, i=$(i), pHasVert[threadid()]=$(pHasVert[threadid()])")
        @assert vState[i][] == threadid() "vState[i][] == threadid()"
        pushDebugMsg!(threadDebugState, "pid=$(Threads.threadid()), PRE-claimAllNbrs")
        #there shouldn't be any other competing threads with current code #RAT
        #RJK 
        claimAllNbrsSucceeded = claimAllNbrs(threadDebugState,pState,a,i,claimStack[threadid()],vState,vLock)
        #assert(claimAllNbrsSucceeded) #there shouldn't be any other competing threads with current code #RAT
        if !claimAllNbrsSucceeded #failure
            myAbandon(threadDebugState,pState[threadid()],vState,vLock,claimStack[threadid()])
            continue  #leave the for loop #TODO: careful if switching to top-level for + inside while
        end
        pushDebugMsg!(threadDebugState, "pid=$(Threads.threadid()), POST-claimAllNbrs")
        # need to check state to switch ELIMINATING or abandon
        pushDebugMsg!(threadDebugState, "pid=$(Threads.threadid()), at PRE-ELIM-STATE-ENTRY")

        selfPsOld = atomic_cas!(pState[threadid()],PS_CLAIMING,PS_ELIMINATING) #0 is for unclaimed
        if selfPsOld != PS_CLAIMING  #someone switched us to assisting
            pushDebugMsg!(threadDebugState, "pid=$(Threads.threadid()), abandoning PRE assert")
            @assert  pState[threadid()][] == PS_ASSISTING "thread assisting"
            pushDebugMsg!(threadDebugState, "pid=$(Threads.threadid()), abandoning POST assist assert")
            myAbandon(threadDebugState, pState[threadid()],vState,vLock,claimStack[threadid()])
            continue  #leave the for loop #TODO: careful if switching to top-level for + inside while
        end
        @assert  pState[threadid()][] == PS_ELIMINATING "thread eliminating"

        pushDebugMsg!(threadDebugState, "pid=$(Threads.threadid()), at POST-ELIM-STATE-ENTRY")

        threadPrintLLcol(a,i)

        len = get_ll_col(a, i, colspace[threadid()])

        pushDebugMsg!(threadDebugState, "pid=$(Threads.threadid()), POST-get_ll_col, len = $(len)")

        #len = compressCol!(a, colspace[threadid()], len, pq, pqLock)  #3hog #TODO #THREADLOCAL #make this threadlocal?
        len = TESTPQcompressCol!(a, colspace[threadid()], len) #TESTPQ #RAT #TODO 

        pushDebugMsg!(threadDebugState, "pid=$(Threads.threadid()), POST-compressCol, len = $(len)")

        csum = zero(Tval) #THREADLOCAL
        for ii in 1:len
            vals[threadid()][ii] = colspace[threadid()][ii].val  #TODO #THREADLOCAL (both vals[threadid()] and colspace[threadid()] need to be made so)
            csum = csum + colspace[threadid()][ii].val  #TODO #THREADLOCAL csum 
            cumspace[threadid()][ii] = csum  #TODO #THREADLOCAL
        end
        wdeg = csum #THREADLOCAL

        pushDebugMsg!(threadDebugState, "pid=$(Threads.threadid()), POST-csum computation ")


        colScale = one(Tval) #THREADLOCAL

        #TODO - we need to switch to manual growing of array?
        local_ldli_row_ptr = atomic_add!(ldli_row_ptr,len) #THREADLOCAL #grab a region of array to use as row for this elim


        pushDebugMsg!(threadDebugState, "pid=$(Threads.threadid()), POST-atomicadd PRE-ldli.colptr[it] alloc")
        #RJK WORKING HERE #TODO
        #MAYBE the way I'm setting 'it' var doesn't make sense on loop-retries?!?!

        pushDebugMsg!(threadDebugState, "pid=$(Threads.threadid()),  it=$(it), local_ldli_row_ptr=$(local_ldli_row_ptr), i=$(i) ")
        pushDebugMsg!(threadDebugState, "pid=$(Threads.threadid()), length(ldli.col)=$(length(ldli.col))")
        ldli.col[it] = i # conversion!  #THREADGLOBAL col
        pushDebugMsg!(threadDebugState, "pid=$(Threads.threadid()), POST col[it]=i")
        ldli.colptr[it] = local_ldli_row_ptr #THREADGLOBAL colptr
        pushDebugMsg!(threadDebugState, "pid=$(Threads.threadid()), POST-colptr ldli.colptr[it] ")
        ldli.collen[it] = len

        for joffset in 1:(len-1)
            pushDebugMsg!(threadDebugState, "pid=$(Threads.threadid()), joffset=$(joffset), local_ldli_row_ptr=$(local_ldli_row_ptr) ")
            ll = colspace[threadid()][joffset] #THREADLOCAL ll  #TODO #THREADGLOBAL CLAIMED
            w = vals[threadid()][joffset] * colScale #THREADLOCAL w  #TODO #THREADGLOBAL CLAIMED?
            j = ll.row #THREADLOCAL j 
            revj = ll.reverse #TODO #THREADGLOBAL CLAIMED?

            f = w/(wdeg) #THREADLOCAL f,w,wdeg

            vals[threadid()][joffset] = zero(Tval) #TODO #THREADLOCAL

            #rked #TODO How does this update for the next probs?
            # think: it should make sense that we zero out things that have been used..
            # AHH, this scheme is not as WO replacement as I thought?
            r = rand(rng) * (csum - cumspace[threadid()][joffset]) + cumspace[threadid()][joffset] #rked, wait, this probably picks a LATER edge? #THREADLOCAL r
            #HIGHTODO rand() IS NOT THREADSAFE!!!
            #TODO rand(rng) is ok?
            koff = searchsortedfirst(cumspace[threadid()],r,one(len),len, o) #TODO is searchsortedfirst threadsafe?!!

            k = colspace[threadid()][koff].row  #TODO #THREADLOCAL koff

            #MUTEX enter mutual exclusion region
            lock(pqLock)
            #approxCholPQInc!(pq, k) #LOWTODO maybe get rid of MUTEX? #TESTPQ
            atomic_fence()
            unlock(pqLock)
            #EXIT mutual exclusion region

            newEdgeVal = f*(one(Tval)-f)*wdeg #THREADLOCAL newEdgeVal

            # fix row k in col j
            revj.row = k   # dense time hog: presumably becaus of cache #TODO #THREADGLOBAL CLAIMED?
            revj.val = newEdgeVal #TODO #THREADGLOBAL CLAIMED?
            revj.reverse = ll #TODO #THREADGLOBAL CLAIMED?

            # fix row j in col k
            khead = a.cols[k] #TODO #THREADGLOBAL CLAIMED?
            a.cols[k] = ll #TODO #THREADGLOBAL CLAIMED?
            ll.next = khead #TODO #THREADGLOBAL CLAIMED?
            ll.reverse = revj #TODO #THREADGLOBAL CLAIMED?
            ll.val = newEdgeVal #TODO #THREADGLOBAL CLAIMED?
            ll.row = j #TODO #THREADGLOBAL CLAIMED?


            colScale = colScale*(one(Tval)-f) #THREADLOCAL colScale
            wdeg = wdeg*(one(Tval)-f)^2 #THREADLOCAL wdeg


            ldli.rowval[local_ldli_row_ptr] = j
            ldli.fval[local_ldli_row_ptr] = f
            #push!(ldli.rowval,j) #TODO #THREADGLOBAL ??
            #push!(ldli.fval, f) #TODO #THREADGLOBAL ??
            #ldli_row_ptr = ldli_row_ptr + one(Tind) #THREADLOCAL #TODO what do we do about this?
            local_ldli_row_ptr += one(Tind) 

            # push!(ops, IJop(i,j,1-f,f))  # another time suck


        end # for joffset

        #LOWTODO broken indentation below        
        ll = colspace[threadid()][len] #TODO #THREADGLOBAL CLAIMED ll?
        w = vals[threadid()][len] * colScale #THREADLOCAL w, colScale #TODO #THREADLOCAL vals[threadid()]
        j = ll.row #THREADLOCAL j #TODO #THREADGLOBAL CLAIMED ll?
        revj = ll.reverse

        pushDebugMsg!(threadDebugState, "pid=$(Threads.threadid()), on elim of i=$(i), final nbr j=$(j)")

        if it < n-1 #OLD COMMENT #rked as we move the it = it+1 update, then this should be it < n-1?
            #MUTEX enter mutual exclusion region
            lock(pqLock)
            atomic_fence()
            #TODO  NEED TO PRINT HERE matlap 13 has "pq.lists[1]=0" to "pq.lists[1]=3" changing here
            pushDebugMsg!(threadDebugState,"pid=$(Threads.threadid()), PRE-DEC-PRINT") #RAT
            #printPQ(pq)
            pushDebugMsg!(threadDebugState,"pid=$(Threads.threadid()), PRE-DEC of j=$(j)") #RAT
            #RJK WORKING HERE threadApproxCholPQDec!(pq, j) #LOWTODO maybe get rid of MUTEX? #TESTPQ
            pushDebugMsg!(threadDebugState,"pid=$(Threads.threadid()), POST-DEC-PRINT") #RAT
            #printPQ(pq)
            atomic_fence()
            unlock(pqLock)
            #EXIT mutual exclusion region
        end

        revj.val = zero(Tval) #TODO #THREADGLOBAL CLAIMED revj?

        ldli.rowval[local_ldli_row_ptr] = j
        ldli.fval[local_ldli_row_ptr] = one(Tval)
        #push!(ldli.rowval,j) #TODO need to fix the way this is accessed?
        #push!(ldli.fval, one(Tval)) #TODO need to fix the way this is accessed?
        #ldli_row_ptr = ldli_row_ptr + one(Tind) #TODO need to fix the way this is accessed?
        #atomic_add!(ldli_row_ptr,one(Tind)) #WE ALREADY INCREMENTED BY FULL AMOUNT

        d[i] = w

        releaseStack(threadDebugState, vState,vLock,claimStack[threadid()])
        vState[i][] == -1 #eliminated #TODO had access bug here?
        pHasVert[threadid()] = 0 # proc is done eliminating
        pState[threadid()][] = PS_ASSISTING #THREADGLOBAL
        #TODO not sure we should really use this move in and out of PS_CLAIMING

        pushDebugMsg!(threadDebugState, "pid=$(Threads.threadid()), at MARK 7, eliminated i=$(i)")
    end #WARNING BAD INDENTATION

    # ccall(:jl_,Void,(Any,), "pid=$(Threads.threadid()) handled for loop iter procId=$(procId), tryWhileCount=$(tryWhileCount)")

end #WARNING BAD INDENTATION

@show atomicCountThreadsContinuing[] #RAT
@show atomicIt[] #RAT


@assert atomicIt[] >= n "err: did not complete enough elims" #is that right?  

#@show it #RAT
#TODO OOPS HERE WE WANT TO END AS USUAL, WHICH I GUESS MEANS
# it = n #TODO change this horrible hack
#ldli.colptr[n] = ldli_row_ptr[] #TODO need to fix the way this is accessed?
#RJK we don't use the last colptr anymore? (and we removed the space for it?!?)

ldli.d = d 

return ldli
end




#=============================================================

The routines that do the solve.

=============================================================#

function LDLsolver(ldli::LDLinv, b::Vector)
    y = copy(b)

    # mu = mean(y)
    # @inbounds for i in eachindex(y)
    #     y[i] = y[i] - mu
    # end #TODO need to explore -- we do this subtraction at start, so presumably we don't need here? want avoid many applications twice in a row
    
    forward!(ldli, y)

    for i in 1:(length(ldli.d))
        #@inbounds for i in 1:(length(ldli.d)) #RAT enable inbounds 
        if ldli.d[i] != 0
            y[i] /= ldli.d[i]
        else
            y[i] = 0
        end    
    end

    backward!(ldli, y)

    mu = mean(y)
    for i in eachindex(y)#RAT enable inbounds 
        #@inbounds for i in eachindex(y)#RAT enable inbounds 
        y[i] = y[i] - mu
    end

    return y
end

function threadLDLsolver(ldli::threadLDLinv, b::Vector)
    y = copy(b)

    # mu = mean(y)
    # @inbounds for i in eachindex(y)
    #     y[i] = y[i] - mu
    # end #TODO need to explore -- we do this subtraction at start, so presumably we don't need here? want avoid many applications twice in a row
    
    threadForward!(ldli, y)

    for i in 1:(length(ldli.d))
        #@inbounds for i in 1:(length(ldli.d)) #RAT enable inbounds 
        if ldli.d[i] != 0
            y[i] /= ldli.d[i]
        else
            y[i] = 0
        end    
    end

    threadBackward!(ldli, y)

    mu = mean(y)
    for i in eachindex(y)#RAT enable inbounds 
        #@inbounds for i in eachindex(y)#RAT enable inbounds 
        y[i] = y[i] - mu
    end

    return y
end


function threadForward!{Tind,Tval}(ldli::threadLDLinv{Tind,Tval}, y::Vector)
    #@show length(ldli.col)
    for ii in 1:length(ldli.col) #RAT enable inbounds 
        #@inbounds for ii in 1:length(ldli.col) #RAT enable inbounds
        ##ccall(:jl_,Void,(Any,), "")
        i = ldli.col[ii]
        ##ccall(:jl_,Void,(Any,), "i=$(i)")


        j0 = ldli.colptr[ii]
        ##ccall(:jl_,Void,(Any,), "j0=$(j0)")
        j1 = j0+ldli.collen[ii]-one(Tind) #TODO FIX 
        ##ccall(:jl_,Void,(Any,), "j1=$(j1)")

        yi = y[i]
        ##ccall(:jl_,Void,(Any,), "yi=$(yi)")
        
        for jj in j0:(j1-1)
            j = ldli.rowval[jj]
            ##ccall(:jl_,Void,(Any,), "j=$(j)")
            y[j] += ldli.fval[jj] * yi
            ##ccall(:jl_,Void,(Any,), "y[j]=$(y[j])")
            yi *= (one(Tval)-ldli.fval[jj])
            ##ccall(:jl_,Void,(Any,), "yi=$(yi)")
        end
        j = ldli.rowval[j1]
        ##ccall(:jl_,Void,(Any,), "j=$(j)")
        y[j] += yi
        ##ccall(:jl_,Void,(Any,), "y[j]=$(y[j])")
        y[i] = yi
        ##ccall(:jl_,Void,(Any,), "yi=$(yi)")
    end
end

function threadBackward!{Tind,Tval}(ldli::threadLDLinv{Tind,Tval}, y::Vector)
    o = one(Tind)
    for ii in length(ldli.col):-1:1
        #@inbounds for ii in length(ldli.col):-1:1 #RAT enable inbounds 
        i = ldli.col[ii]

        j0 = ldli.colptr[ii]
        j1 = j0+ldli.collen[ii]-o #TODO TEST 

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

# function forward!{Tind,Tval}(ldli::LDLinv{Tind,Tval}, y::Vector)
#     #@show length(ldli.col)
#     for ii in 1:length(ldli.col) #RAT enable inbounds 
#         #@inbounds for ii in 1:length(ldli.col) #RAT enable inbounds
#         #ccall(:jl_,Void,(Any,), "")
#         i = ldli.col[ii]
#         #ccall(:jl_,Void,(Any,), "i=$(i)")


#         j0 = ldli.colptr[ii]
#         #ccall(:jl_,Void,(Any,), "j0=$(j0)")
#         j1 = ldli.colptr[ii+1]-one(Tind) #TODO FIX 
#         #ccall(:jl_,Void,(Any,), "j1=$(j1)")

#         yi = y[i]
#         #ccall(:jl_,Void,(Any,), "yi=$(yi)")

#         for jj in j0:(j1-1)
#             j = ldli.rowval[jj]
#             #ccall(:jl_,Void,(Any,), "j=$(j)")
#             y[j] += ldli.fval[jj] * yi
#             #ccall(:jl_,Void,(Any,), "y[j]=$(y[j])")
#             yi *= (one(Tval)-ldli.fval[jj])
#             #ccall(:jl_,Void,(Any,), "yi=$(yi)")
#         end
#         j = ldli.rowval[j1]
#         #ccall(:jl_,Void,(Any,), "j=$(j)")
#         y[j] += yi
#         #ccall(:jl_,Void,(Any,), "y[j]=$(y[j])")
#         y[i] = yi
#         #ccall(:jl_,Void,(Any,), "yi=$(yi)")
#     end
# end

function forward!{Tind,Tval}(ldli::LDLinv{Tind,Tval}, y::Vector)
    #@show length(ldli.col)
    for ii in 1:length(ldli.col) #RAT enable inbounds 
        #@inbounds for ii in 1:length(ldli.col) #RAT enable inbounds
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

function backward!{Tind,Tval}(ldli::LDLinv{Tind,Tval}, y::Vector)
    o = one(Tind)
    for ii in length(ldli.col):-1:1
        #@inbounds for ii in length(ldli.col):-1:1 #RAT enable inbounds 
        i = ldli.col[ii]

        j0 = ldli.colptr[ii]
        j1 = ldli.colptr[ii+1]-o #TODO FIX 

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
print#y[:,i] = yi
end
end

=#


"""
solver = approxCholLap(a; tol::Real=1e-6, maxits=1000, maxtime=Inf, verbose=false, pcgIts=Int[], params=ApproxCholParams())

A heuristic by Daniel Spielman inspired by the linear system solver in https://arxiv.org/abs/1605.02353 by Rasmus Kyng and Sushant Sachdeva.  Whereas that paper eliminates vertices one at a time, this eliminates edges one at a time.  It is probably possible to analyze it.
The `ApproxCholParams` let you choose one of three orderings to perform the elimination.

* ApproxCholParams(:given) - in the order given.
This is the fastest for construction the preconditioner, but the slowest solve.
* ApproxCholParams(:deg) - always eliminate the node of lowest degree.
This is the slowest build, but the fastest solve.
* ApproxCholParams(:wdeg) - go by a perturbed order of wted degree.
This is the sweet spot in between.
"""
function approxCholLap{Tv,Ti}(a::SparseMatrixCSC{Tv,Ti};
  tol::Real=1e-6,
  maxits=1000,
  maxtime=Inf,
  verbose=false,
  pcgIts=Int[],
  params=ApproxCholParams())

return lapWrapComponents(approxCholLap1, a,
   verbose=verbose,
   tol=tol,
   maxits=maxits,
   maxtime=maxtime,
   pcgIts=pcgIts,
   params=params)


end

function approxCholLapGreedy(a::SparseMatrixCSC;
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

llmat = LLmatp(a)
ldli = approxChol(llmat)
F(b) = LDLsolver(ldli, b)

if verbose
    println("Using greedy degree ordering. Factorization time: ", time()-t1)
    println("Ratio of operator edges to original edges: ", 2 * length(ldli.fval) / nnz(a))
end

if verbose
    println("ratio of max to min diagonal of laplacian : ", maximum(diag(la))/minimum(diag(la)))
end


f(b;tol=tol_,maxits=maxits_, maxtime=maxtime_, verbose=verbose_, pcgIts=pcgIts_) = pcg(la, b-mean(b), F, tol=tol, maxits=maxits, maxtime=maxtime, pcgIts=pcgIts, verbose=verbose, stag_test = params.stag_test)

end

function approxCholLapGiven(a::SparseMatrixCSC;
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
ldli = approxChol(llmat)
F(b) = LDLsolver(ldli, b)

if verbose
    println("Using given ordering. Factorization time: ", time()-t1)
    println("Ratio of operator edges to original edges: ", 2 * length(ldli.fval) / nnz(a))
end

if verbose
    println("ratio of max to min diagonal of laplacian : ", maximum(diag(la))/minimum(diag(la)))
end


f(b;tol=tol_,maxits=maxits_, maxtime=maxtime_, verbose=verbose_, pcgIts=pcgIts_) = pcg(la, b-mean(b), F, tol=tol, maxits=maxits, maxtime=maxtime, pcgIts=pcgIts, verbose=verbose, stag_test = params.stag_test)

end

function approxCholLapWdeg(a::SparseMatrixCSC;
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

v = vec(sum(a,1))
v = v .* (1+ rand(length(v)))
p = sortperm(v)

llmat = LLMatOrd(a,p)
ldli = approxChol(llmat)

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



function approxCholLap1{Tv,Ti}(a::SparseMatrixCSC{Tv,Ti};
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


if params.order == :deg

    return approxCholLapGreedy(a,
     verbose=verbose,
     tol=tol,
     maxits=maxits,
     maxtime=maxtime,
     pcgIts=pcgIts,
     params=params)


elseif params.order == :wdeg

    return approxCholLapWdeg(a,
       verbose=verbose,
       tol=tol,
       maxits=maxits,
       maxtime=maxtime,
       pcgIts=pcgIts,
       params=params)


else
    return approxCholLapGiven(a,
      verbose=verbose,
      tol=tol,
      maxits=maxits,
      maxtime=maxtime,
      pcgIts=pcgIts,
      params=params)


end

end

# """
#     solver = approxCholSddm(sddm; tol::Real=1e-6, maxits=1000, maxtime=Inf, verbose=false, pcgIts=Int[], params=ApproxCholParams())

# Solves sddm systems by wrapping approxCholLap.
# Not yet optimized directly for sddm.
# """
# approxCholSddm = sddmWrapLap(approxCholLap)




    #===============================

    Checking the condition number

    =================================#

    """
    cn = condNumber(a, ldli; verbose=false)

    Given an adjacency matrix a and an ldli computed by approxChol,
    this computes the condition number.
    """
    function condNumber(a, ldli; verbose=false)
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
        upper = eigs(gOp;nev=1,which=:LM,tol=1e-2)[1][1]

        g2(b) = upper*b - g(b)
        g2Op = SqLinOp(true,1.0,size(a,1),g2)
        lower = upper - eigs(g2Op;nev=2,which=:LM,tol=1e-2)[1][2]

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
        This variation of approxChol creates a cholesky factor to do the elimination.
        It has not yet been optimized, and does not yet make the cholesky factor lower triangular
        """
        function approxCholLapChol{Tv,Ti}(a::SparseMatrixCSC{Tv,Ti}; tol::Real=1e-6, maxits=1000, maxtime=Inf, verbose=false, pcgIts=Int[])

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


            f(b;tol=tol_,maxits=maxits_, maxtime=maxtime_, verbose=verbose_, pcgIts=pcgIts_) = pcg(la, b-mean(b), F, tol=tol, maxits=maxits, maxtime=maxtime, pcgIts=pcgIts, verbose=verbose)

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

function printPQ{Tind}(pq::ThreadApproxCholPQ{Tind})
    ccall(:jl_,Void,(Any,), "pid=$(Threads.threadid()) printing PQ")
    ccall(:jl_,Void,(Any,), "pid=$(Threads.threadid()) pq.nitems=$(pq.nitems)")
    ccall(:jl_,Void,(Any,), "pid=$(Threads.threadid()) pq.nitemsAtomic=$(pq.nitemsAtomic[])")
    for k in eachindex(pq.elems)
        ccall(:jl_,Void,(Any,), "pid=$(Threads.threadid()) pq.elems[$(k)]=$(pq.elems[k])")
    end
    for l in eachindex(pq.lists)
        ccall(:jl_,Void,(Any,), "pid=$(Threads.threadid()) pq.lists[$(l)]=$(pq.lists[l])")
    end
    ccall(:jl_,Void,(Any,), "pid=$(Threads.threadid()) pq.minlist=$(pq.minlist)")
end

function printDebug{Tind}(tds::ThreadDebugState{Tind})
    ccall(:jl_,Void,(Any,), "pid=$(Threads.threadid()) DEBUG PRINT")
    #tds. #HIGHTODO this must be there by mistake -- or did I stop in the middle of something?
    atomic_fence()
    printPQ(tds.pq)
    for p = 1:nthreads()
        ccall(:jl_,Void,(Any,), "pid=$(Threads.threadid()) DEBUG PRINT p=$(p), pState[p]=$(tds.pState[p]), pHasVert[p]=$(tds.pHasVert[p])")
    end
    ccall(:jl_,Void,(Any,), "pid=$(Threads.threadid()) DEBUG PRINT MESSAGE QUEUES")
    for p = 1:nthreads()
        for m = 1:length(tds.pDebugMsgQueue[p])
            ccall(:jl_,Void,(Any,),tds.pDebugMsgQueue[p][m])
            
        end
    end
end

function pushDebugMsg!{Tind}(tds::ThreadDebugState{Tind},msg::String)
    atomic_fence() #TODO this might not work well
    push!(tds.pDebugMsgQueue[threadid()],msg)
end


# #TODO are these of acceptable types?! probably not
# const PS_ASSISTING = convert(Int64,0)
# const PS_CLAIMING = convert(Int64,1)
# const PS_ELIMINATING = convert(Int64,2)


function ApproxCholPQ{Tind}(a::Vector{Tind})
    n = length(a)
    elems = Array{ApproxCholPQElem{Tind}}(n)
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

    return ApproxCholPQ(elems, lists, minlist, n, Atomic{Int64}(n), n)
end

function ThreadApproxCholPQ{Tind}(a::Vector{Tind})

    n = length(a)
    elems = Array{ApproxCholPQElem{Tind}}(n)
    lists = zeros(Tind, 2*n+1)
    inQueue = ones(Bool,n)
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

    return ThreadApproxCholPQ(elems, lists, minlist, inQueue, n, Atomic{Int64}(n), n)
end

function approxCholPQPop!{Tind}(pq::ApproxCholPQ{Tind})
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
    prevVal = Threads.atomic_sub!(pq.nitemsAtomic,1)
    @assert pq.nitems == prevVal-1 "item counts out of sync" 

    #pq.lists[9] = rand(Int64) #RAT
    
    return i
end

function threadApproxCholPQPop!{Tind}(pq::ThreadApproxCholPQ{Tind})
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
    prevVal = Threads.atomic_sub!(pq.nitemsAtomic,1)
    @assert pq.nitems == prevVal-1 "item counts out of sync" 

    #pq.lists[9] = rand(Int64) #RAT

    pq.inQueue[i] = false
    
    return i
end

function approxCholPQMove!{Tind}(pq::ApproxCholPQ{Tind}, i, newkey, oldlist, newlist)

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

    return Void
end

function threadApproxCholPQMove!{Tind}(pq::ThreadApproxCholPQ{Tind}, i, newkey, oldlist, newlist)

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

    return Void
end

"""
Decrement the key of element i
This could crash if i exceeds the maxkey
"""
function approxCholPQDec!{Tind}(pq::ApproxCholPQ{Tind}, i)

    oldlist = keyMap(pq.elems[i].key, pq.n)
    newlist = keyMap(pq.elems[i].key - one(Tind), pq.n)

    if newlist != oldlist

        approxCholPQMove!(pq, i, pq.elems[i].key - one(Tind), oldlist, newlist)

        if newlist < pq.minlist
            pq.minlist = newlist
        end

    else
        pq.elems[i] = ApproxCholPQElem{Tind}(pq.elems[i].prev, pq.elems[i].next, pq.elems[i].key - one(Tind))
    end


    return Void
end

function threadApproxCholPQDec!{Tind}(pq::ThreadApproxCholPQ{Tind}, i)
    if !pq.inQueue[i]
        return Void
    end

    oldlist = keyMap(pq.elems[i].key, pq.n)
    newlist = keyMap(pq.elems[i].key - one(Tind), pq.n)

    if newlist != oldlist

        threadApproxCholPQMove!(pq, i, pq.elems[i].key - one(Tind), oldlist, newlist)

        if newlist < pq.minlist
            pq.minlist = newlist
        end

    else
        pq.elems[i] = ApproxCholPQElem{Tind}(pq.elems[i].prev, pq.elems[i].next, pq.elems[i].key - one(Tind))
    end


    return Void
end

"""
Increment the key of element i
This could crash if i exceeds the maxkey
"""
function approxCholPQInc!{Tind}(pq::ApproxCholPQ{Tind}, i)

    oldlist = keyMap(pq.elems[i].key, pq.n)
    newlist = keyMap(pq.elems[i].key + one(Tind), pq.n)

    if newlist != oldlist

        approxCholPQMove!(pq, i, pq.elems[i].key + one(Tind), oldlist, newlist)

    else
        pq.elems[i] = ApproxCholPQElem{Tind}(pq.elems[i].prev, pq.elems[i].next, pq.elems[i].key + one(Tind))
    end

    return Void
end

function threadApproxCholPQInc!{Tind}(pq::ThreadApproxCholPQ{Tind}, i)
    if !pq.inQueue[i]
        return Void
    end
    
    oldlist = keyMap(pq.elems[i].key, pq.n)
    newlist = keyMap(pq.elems[i].key + one(Tind), pq.n)

    if newlist != oldlist

        threadApproxCholPQMove!(pq, i, pq.elems[i].key + one(Tind), oldlist, newlist)

    else
        pq.elems[i] = ApproxCholPQElem{Tind}(pq.elems[i].prev, pq.elems[i].next, pq.elems[i].key + one(Tind))
    end

    return Void
end
