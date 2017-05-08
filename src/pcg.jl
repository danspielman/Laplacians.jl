#=

Implementations of cg and pcg that use BLAS when they can.

Started by Dan Spielman.
Contributors:

=#


"""
    x = cg(mat, b; tol, maxits, maxtime, verbose, pcgIts)

solves a symmetric linear system `mat x = b`.
# Arguments
* `tol` is set to 1e-6 by default,
* `maxits` defaults to Inf
* `maxtime` defaults to Inf.  It measures seconds.
* `verbose` defaults to false
* `pcgIts` is an array for returning the number of pcgIterations.  Default is length 0, in which case nothing is returned.
"""
function cg end

"""
    x = pcg(mat, b, pre; tol, maxits, maxtime, verbose, pcgIts)` 

solves a symmetric linear system using preconditioner `pre`.
# Arguments
* `pre` can be a function or a matrix.  If a matrix, a function to solve it is created with cholFact.
* `tol` is set to 1e-6 by default,
* `maxits` defaults to Inf
* `maxtime` defaults to Inf.  It measures seconds.
* `verbose` defaults to false
* `pcgIts` is an array for returning the number of pcgIterations.  Default is length 0, in which case nothing is returned.
"""
function pcg end


"""
    x = cgSolver(mat; tol, maxits, maxtime, verbose, pcgIts)

creates a solver for a PSD system `mat`.
The parameters are as described in cg.
"""
function cgSolver end


"""
    x = pcgSolver(mat, pre; tol, maxits, maxtime, verbose, pcgIts)

creates a solver for a PSD system using preconditioner `pre`.
The parameters are as described in pcg.
"""
function pcgSolver end

function cgSolver(mat; tol::Real=1e-6, maxits=Inf, maxtime=Inf, verbose=false, pcgIts=Int[])

    tol_=tol
    maxits_=maxits
    maxtime_=maxtime
    verbose_=verbose
    pcgIts_=pcgIts


    f(b; tol=tol_, maxits=maxits_, maxtime=maxtime_, verbose=verbose_, pcgIts=pcgIts_) =
       cg(mat, b, tol=tol, maxits=maxits, maxtime=maxtime, verbose=verbose, pcgIts=pcgIts)
end
    

"""
    x = cgLapSolver(A::AbstractMatrix; tol::Real=1e-6, maxits=Inf, maxtime=Inf, verbose=false, pcgIts=Int[])

Create a solver that uses cg to solve Laplacian systems in the laplacian of A. 
This just exists to satisfy our interface.
It does nothing more than create the Laplacian and call cg on each connected component.
"""
function cgLapSolver(a::SparseMatrixCSC; tol::Real=1e-6, maxits=Inf, maxtime=Inf, verbose=false, pcgIts=Int[])

    return lapWrapComponents(cgLapSolver1, a, verbose=verbose, tol=tol, maxits=maxits, maxtime=maxtime, pcgIts=pcgIts)


end


function cgLapSolver1(A::AbstractMatrix; tol::Real=1e-6, maxits=Inf, maxtime=Inf, verbose=false, pcgIts=Int[])

    la = forceLap(A)
    
    tol_=tol
    maxits_=maxits
    maxtime_=maxtime
    verbose_=verbose
    pcgIts_=pcgIts

    f(b; tol=tol_, maxits=maxits_, maxtime=maxtime_, verbose=verbose_, pcgIts=pcgIts_) =
        cg(la, b-mean(b), tol=tol, maxits=maxits, maxtime=maxtime, verbose=verbose, pcgIts=pcgIts)

end



function pcg(mat, b, pre::Union{AbstractArray,Matrix}; tol::Real=1e-6, maxits=Inf, maxtime=Inf, verbose=false, pcgIts=Int[])
    fact = cholfact(pre)
    F = x->(fact \ x)
    pcg(mat, b, F; tol=tol, maxits=maxits, maxtime=maxtime, verbose=verbose, pcgIts=pcgIts)
end


function pcgSolver(mat, pre; tol::Real=1e-6, maxits=Inf, maxtime=Inf, verbose=false, pcgIts=Int[])
    tol_=tol
    maxits_=maxits
    maxtime_=maxtime
    verbose_=verbose
    pcgIts_=pcgIts

    f(b; tol=tol_, maxits=maxits_, maxtime=maxtime_, verbose=verbose_, pcgIts=pcgIts_) =
        pcg(mat, b, pre, tol=tol, maxits=maxits, maxtime=maxtime, verbose=verbose, pcgIts=pcgIts)
end


"""
    x = pcgLapSolver(A, B; tol::Real=1e-6, maxits=Inf, maxtime=Inf, verbose=false, pcgIts=Int[])

Create a solver that uses pcg to solve Laplacian systems in `A`
Specialized for the case when the preconditioner the Laplacian matrix of `B`.
It solves the preconditioner by Cholesky Factorization.
"""
function pcgLapSolver(A::AbstractMatrix, B::AbstractMatrix; tol::Real=1e-6, maxits=Inf, maxtime=Inf, verbose=false, pcgIts=Int[])
    fact = cholLap(B)
    
    tol_=tol
    maxits_=maxits
    maxtime_=maxtime
    verbose_=verbose
    pcgIts_=pcgIts

    la = forceLap(A)

    f(b; tol=tol_, maxits=maxits_, maxtime=maxtime_, verbose=verbose_, pcgIts=pcgIts_) =
        pcg(la, b, fact, tol=tol, maxits=maxits, maxtime=maxtime, verbose=verbose, pcgIts=pcgIts)

end

function cg{Tval}(mat, b::Array{Tval,1};
        tol::Real=1e-6, maxits=Inf, maxtime=Inf, verbose=false, pcgIts=Int[])

    local al::Tval

    n = size(mat,2)

    nb = norm(b)

    # If input vector is zero, quit
    if nb == 0
      return bestx
    end

    x = zeros(Tval,n)
    bestx = zeros(Tval,n)
    bestnr = one(Tval)

    r = copy(b)

    p = copy(r)

    rho = norm(r)^2

    t1 = time()

    itcnt = 0
    while itcnt < maxits
        itcnt = itcnt+1

        q = mat*p

        al = rho/dot(p, q)

        axpy2!(al,p,x)
        # x = x + al * p

        axpy2!(-al,q,r)
        #r .= r .- al.*q

        nr = norm(r)/nb
        if nr < bestnr
          bestnr = nr
          @inbounds @simd for i in 1:n
            bestx[i] = x[i]
          end
        end
        if nr < tol #Converged?
            break
        end

        # here is the top of the code in numerical templates

        oldrho = rho
        rho = norm(r)^2
        if (rho < eps(Tval) || isinf(rho))
            println("CG stopped for rho: ", rho)
          break
        end

        # the following would have higher accuracy
        #       rho = sum(r.^2)

        beta = rho/oldrho
        if (beta< eps(Tval) || isinf(beta))
            println("CG stopped for beta: ", beta)

          break
        end

        bzbeta!(beta,p,r)
        # p = z + beta*p

        if (time() - t1) > maxtime
            if verbose
                println("CG stopped at maxtime.")
            end
            break
        end

    end

    if verbose
        println("CG stopped after: ", round((time() - t1),3), " seconds and ", itcnt, " iterations with relative error ", (norm(r)/norm(b)), ".")
    end

    if length(pcgIts) > 0
        pcgIts[1] = itcnt
    end


    return bestx
end

function pcg{Tval}(mat, b::Array{Tval,1}, pre::Function;
        tol::Real=1e-6, maxits=Inf, maxtime=Inf, verbose=false, pcgIts=Int[])

    local al::Tval

    n = size(mat,2)

    nb = norm(b)

    # If input vector is zero, quit
    if nb == 0
      return bestx
    end

    x = zeros(Tval,n)
    bestx = zeros(Tval,n)
    bestnr = one(Tval)

    r = copy(b)
    z = pre(r)
    p = copy(z)

    rho = dot(r, z)

    t1 = time()

    itcnt = 0
    while itcnt < maxits
        itcnt = itcnt+1

        q = mat*p

        al = rho/dot(p, q)

        axpy2!(al,p,x)
        # x = x + al * p
        #=
        @inbounds @simd for i in 1:n
            x[i] += al*p[i]
        end
        =#
        #axpy

        axpy2!(-al,q,r)
        #r .= r .- al.*q
        #=
        @inbounds @simd for i in 1:n
            r[i] -= al*q[i]
        end
        =#

        nr = norm(r)/nb
        if nr < bestnr
          bestnr = nr
          @inbounds @simd for i in 1:n
            bestx[i] = x[i]
          end
        end
        if nr < tol #Converged?
            break
        end

        # here is the top of the code in numerical templates

        z = pre(r)

        oldrho = rho
        rho = dot(z, r)
        if (rho < eps(Tval) || isinf(rho))
          break
        end

        # the following would have higher accuracy
        #       rho = sum(r.^2)

        beta = rho/oldrho
        if (beta< eps(Tval) || isinf(beta))
          break
        end

        bzbeta!(beta,p,z)
        #=
        # p = z + beta*p
        @inbounds @simd for i in 1:n
            p[i] = z[i] + beta*p[i]
        end
        =#

        if (time() - t1) > maxtime
            if verbose
                println("PCG New stopped at maxtime.")
            end
            break
        end

    end

    if verbose
        println("PCG stopped after: ", round((time() - t1),3), " seconds and ", itcnt, " iterations with relative error ", (norm(r)/norm(b)), ".")
    end

    if length(pcgIts) > 0
        pcgIts[1] = itcnt
    end


    return bestx
end


function axpy2!(al,p::Array,x::Array)
  n = length(x)
  @inbounds @simd for i in 1:n
      x[i] = x[i] + al*p[i]
  end
end

# p[i] = z[i] + beta*p[i]
function bzbeta!(beta,p::Array,z::Array)
  n = length(p)
  @inbounds @simd for i in 1:n
      p[i] = z[i] + beta*p[i]
  end
end






#==========================================================
   Code we no longer use
===========================================================#
   

# uses BLAS.  As fast as Matlab's pcg.
# set to use similar paramaters
function cgBLAS{Tval}(mat, b::Array{Tval,1};
        tol::Real=1e-6, maxits=Inf, maxtime=Inf, verbose=false, pcgIts=Int[])

    n = size(mat,2)
    
    x = zeros(Tval,n)
    
    tol = tol * norm(b)
    r = copy(b) 
    p = copy(r)

    # the following would have higher accuracy
    #       rho = sum(r.^2)
    rho = norm(r)^2

    t1 = time()

    itcnt = 0
    while itcnt < maxits
        itcnt = itcnt+1
    
        q = mat*p

        al = rho/dot(p, q)

        BLAS.axpy!(al,p,x)  # x = x + al * p

        BLAS.axpy!(-al,q,r)  # r -= al*q

        if norm(r) < tol #Converged?
            break
        end

        oldrho = rho
        rho = norm(r)^2

        # the following would have higher accuracy
        #       rho = sum(r.^2)
        
        beta = rho/oldrho

        # p = r + beta*p
        BLAS.scal!(n,beta,p,1) # p *= beta
        BLAS.axpy!(1.0,r,p) # p += r

       
        if (time() - t1) > maxtime
            if verbose
                println("CG BLAS stopped at maxtime.")
            end
            break
        end
      end

    if verbose
        println("CG BLAS stopped after: ", round((time() - t1),3), " seconds and ", itcnt, " iterations with relative error ", (norm(r)/norm(b)), ".")

    end

    if length(pcgIts) > 0
        pcgIts[1] = itcnt
    end
    
    return x
end

# For when you can't use BLAS.
# set to use similar paramaters similar to MATLAB
function cgSlow{Tval}(mat, b::Array{Tval,1};
        tol::Real=1e-6, maxits=Inf, maxtime=Inf, verbose=false, pcgIts=Int[])

    n = size(mat,2)
    
    x = zeros(Tval,n)
    
    tol = tol * norm(b)
    r = copy(b) 
    p = copy(r)

    # the following would have higher accuracy
    #       rho = sum(r.^2)
    rho = norm(r)^2

    t1 = time()
    
    itcnt = 0
    while itcnt < maxits
        itcnt = itcnt+1
    
        q = mat*p

        al = rho/dot(p, q)

        # x = x + al * p
        for i in 1:n
            x[i] = x[i] + al*p[i]
        end

        # r = r - al*q
        for i in 1:n
            r[i] = r[i] - al*q[i]
        end

        if norm(r) < tol 
            break
        end

        oldrho = rho
        rho = norm(r)^2

        # the following would have higher accuracy
        #       rho = sum(r.^2)
        
        beta = rho/oldrho


        # p = r + beta*p
        for i in 1:n
            p[i] = r[i] + beta*p[i]
        end

        if (time() - t1) > maxtime
            if verbose
                println("CG Slow stopped at maxtime.")
            end
            break
        end
      end

    if verbose
        println("CG Slow stopped after: ", round((time() - t1),3), " seconds and ", itcnt, " iterations with relative error ", (norm(r)/norm(b)), ".")
    end

    if length(pcgIts) > 0
        pcgIts[1] = itcnt
    end

    return x
end

#====================================
   PCG Code
=====================================#


# uses BLAS.  As fast as Matlab's pcg.
# set to use similar parameters "
function pcgBLAS{Tval}(mat, b::Array{Tval,1}, pre;
        tol::Real=1e-6, maxits=Inf, maxtime=Inf, verbose=false, pcgIts=Int[])

    n = size(mat,2)
    
    x = zeros(Tval,n)
    
    tol = tol * norm(b)
    r = copy(b) 
    z = pre(r)
    p = copy(z)

    rho = dot(r, z) # BLAS.dot does not seem faster

    t1 = time()

    itcnt = 0
    while itcnt < maxits
        itcnt = itcnt+1
        
        q = mat*p

        if dot(p,q) == 0
            break
        end

        al = rho/dot(p, q)

        BLAS.axpy!(al,p,x)  # x = x + al * p

        BLAS.axpy!(-al,q,r)  # r -= al*q

        if norm(r) < tol #Converged?
            break
        end

        # here is the top of the code in numerical templates

        z = pre(r)

        oldrho = rho
        rho = dot(z, r)

        # the following would have higher accuracy
        #       rho = sum(r.^2)
        
        beta = rho/oldrho



        p = z + beta*p
        # BLAS.scal!(n,beta,p,1) # p *= beta
        # BLAS.axpy!(1.0,z,p) # p += z

        if (time() - t1) > maxtime
            if verbose
                println("PCG BLAS stopped at maxtime.")
            end
            break
        end

       
      end

    if verbose
        println("PCG BLAS stopped after: ", round((time() - t1),3), " seconds and ", itcnt, " iterations with relative error ", (norm(r)/norm(b)), ".")
    end

    if length(pcgIts) > 0
        pcgIts[1] = itcnt
    end


    return x
end

function pcgSlow{Tval}(mat, b::Array{Tval,1}, pre;
        tol::Real=1e-6, maxits=Inf, maxtime=Inf, verbose=false, pcgIts=Int[])


    n = size(mat,2)
    
    x = zeros(Tval,n)
    
    tol = tol * norm(b)
    r = copy(b) 
    z = pre(r)
    p = copy(z)

    rho = dot(r, z) 

    t1 = time()

    itcnt = 0
    while itcnt < maxits
        itcnt = itcnt+1
    
        q = mat*p

        al = rho/dot(p, q)

        # x = x + al * p
        for i in 1:n
            x[i] = x[i] + al*p[i]
        end

        # r = r - al*q
        for i in 1:n
            r[i] = r[i] - al*q[i]
        end

        if norm(r) < tol #Converged?
            break
        end

        # here is the top of the code in numerical templates

        z = pre(r)

        oldrho = rho
        rho = dot(z, r)

        # the following would have higher accuracy
        #       rho = sum(r.^2)
        
        beta = rho/oldrho

        # p = z + beta*p
        for i in 1:n
            p[i] = z[i] + beta*p[i]
        end

        if (time() - t1) > maxtime
            if verbose
                println("PCG Slow stopped at maxtime.")
            end
            break
        end
       
    end

    if verbose
        println("PCG Slow stopped after: ", round((time() - t1),3), " seconds and ", itcnt, " iterations with relative error ", (norm(r)/norm(b)), ".")
    end

    if length(pcgIts) > 0
        pcgIts[1] = itcnt
    end


    return x
end
