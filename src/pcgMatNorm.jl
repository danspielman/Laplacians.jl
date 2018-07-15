 #=

Implementations of pcg that use BLAS when they can.
Computes error in matrix norm, if correct left hand side is present.

Started by Dan Spielman.
Contributors: Serban Stan

=#

"""
`pcg(mat, b, pre; tol, maxits, maxtime, verbose)` solves a symmetric linear system using preconditioner `pre`.
`pre` can be a function or a matrix.  If a matrix, a function to solve it is created with cholFact.
`tol` is set to 1e-6 by default,
`maxits` defaults to Inf
`maxtime` defaults to Inf.  It measures seconds.
`verbose` defaults to false
`lhs` defaults to the all zeros vector
"""

function pcgMatNorm(mat, b::Array{Tv,1}, pre::Function, lhs::Array{Tv,1}; tol::Real=1e-6, maxits=Inf, maxtime=Inf, verbose=false) where Tv
    return pcgBLASMatNorm(mat, b, pre, lhs, tol=tol, maxits=maxits, maxtime=maxtime, verbose=verbose)
end


#====================================
   PCG Code
=====================================#


# uses BLAS.  As fast as Matlab's pcg.
# set to use similar paramaters
function pcgBLASMatNorm(mat, b::Array{Tval,1}, pre, lhs::Array{Tval,1}; tol::Real=1e-6, maxits=Inf, maxtime=Inf, verbose=false) where Tval

    n = size(mat,2)
    
    x = zeros(Tval,n)
    
    r = copy(b) 
    z = pre(r)
    p = copy(z)

    rho = dot(r, z) # BLAS.dot does not seem faster

    t1 = time()

    itcnt = 0
    while itcnt < maxits
        itcnt = itcnt+1
        
        q = mat*p

        al = rho/dot(p, q)

        BLAS.axpy!(al,p,x)  # x = x + al * p
        BLAS.axpy!(-al,q,r)  # r -= al*q

        # This will take some time. Can we do something faster?
        if sqrt((lhs - x)' * mat * (lhs - x))[1] < tol 
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
                println("PCG stopped at maxtime.")
            end
            break
        end

      end

    if verbose
        println("PCG stopped after: ", itcnt, " iterations with relative error ", sqrt((lhs - x)' * mat * (lhs - x))[1], ".")
    end

    return x
end
