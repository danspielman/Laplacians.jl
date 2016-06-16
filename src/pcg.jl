#=

Implementations of cg and pcg that use BLAS when they can.

Started by Dan Spielman.
Contributors:

=#


"""
`cg(mat, b; tol, maxits, maxtime, verbose)` solves a symmetric linear system.
`tol` is set to 1e-6 by default,
`maxits` defaults to Inf
`maxtime` defaults to Inf.  It measures seconds.
`verbose` defaults to false
"""
function cg end

"""
`pcg(mat, b, pre; tol, maxits, maxtime, verbose)` solves a symmetric linear system using preconditioner `pre`.
`pre` can be a function or a matrix.  If a matrix, a function to solve it is created with cholFact.
`tol` is set to 1e-6 by default,
`maxits` defaults to Inf
`maxtime` defaults to Inf.  It measures seconds.
`verbose` defaults to false
"""
function pcg end

function cg(mat, b::Array{Float64,1}; tol::Real=1e-6, maxits=Inf, maxtime=Inf, verbose=false)
    cgBLAS(mat, b, tol=tol, maxits=maxits, maxtime=maxtime, verbose=verbose)
end

function cg(mat, b::Array{Float32,1}; tol::Real=1e-6, maxits=Inf, maxtime=Inf, verbose=false)
    cgBLAS(mat, b, tol=tol, maxits=maxits, maxtime=maxtime, verbose=verbose)
end

function cg(mat, b; tol::Real=1e-6, maxits=Inf, maxtime=Inf, verbose=true)
    cgSlow(mat, b, tol=tol, maxits=maxits, maxtime=maxtime, verbose=verbose)
end

"""Create a solver that uses cg to solve systems in mat.  Fix the default parameters of the solver as given"""
function cgSolver(mat; tol::Real=1e-6, maxits=Inf, maxtime=Inf, verbose=false)
    f(b; tol::Real=1e-6, maxits=maxits, maxtime=maxtime, verbose=verbose) = cg(mat, b, tol=tol, maxits=maxits, maxtime=maxtime, verbose=verbose)
end
    


function pcg(mat, b, pre::Union{AbstractArray,Matrix}; tol::Real=1e-6, maxits=Inf, maxtime=Inf, verbose=false)
    fact = cholfact(pre)
    F = x->(fact \ x)
    pcg(mat, b, F; tol=tol, maxits=maxits, maxtime=maxtime, verbose=verbose)
end


function pcg(mat, b::Array{Float64,1}, pre::Function; tol::Real=1e-6, maxits=Inf, maxtime=Inf, verbose=false)
    pcgBLAS(mat, b, pre, tol=tol, maxits=maxits, maxtime=maxtime, verbose=verbose)
end

function pcg(mat, b::Array{Float32,1}, pre::Function; tol::Real=1e-6, maxits=Inf, maxtime=Inf, verbose=false)
    pcgBLAS(mat, b, pre, tol=tol, maxits=maxits, maxtime=maxtime, verbose=verbose)
end



function pcg(mat, b, pre::Function; tol::Real=1e-6, maxits=Inf, maxtime=Inf, verbose=false)
    pcgSlow(mat, b, pre, tol=tol, maxits=maxits, maxtime=maxtime, verbose=verbose)
end


"""Create a solver that uses pcg to solve systems in mat.  Fix the default parameters of the solver as given"""
function pcgSolver(mat, pre; tol::Real=1e-6, maxits=Inf, maxtime=Inf, verbose=false)
    f(b; tol::Real=1e-6, maxits=maxits, maxtime=maxtime, verbose=verbose) = pcg(mat, b, pre, tol=tol, maxits=maxits, maxtime=maxtime, verbose=verbose)

end

"""Create a solver that uses cg to solve Laplacian systems in mat. Specialized for the case when pre is a Laplacian matrix.  Fix the default parameters of the solver as given"""
function pcgLapSolver(mat, pre; tol::Real=1e-6, maxits=Inf, maxtime=Inf, verbose=false)
    fact = lapChol(pre)
    f(b; tol::Real=1e-6, maxits=maxits, maxtime=maxtime, verbose=verbose) = pcg(mat, b, fact, tol=tol, maxits=maxits, maxtime=maxtime, verbose=verbose)

end





# uses BLAS.  As fast as Matlab's pcg.
# set to use similar paramaters
function cgBLAS{Tval}(mat, b::Array{Tval,1};
        tol::Real=1e-6, maxits=Inf, maxtime=Inf, verbose=false)

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
                println("CG stopped at maxtime.")
            end
            break
        end
      end

    if verbose
        println("CG stopped after: ", itcnt, " iterations with relative error ", (norm(r)/norm(b)), ".")
    end

    return x
end

# For when you can't use BLAS.
# set to use similar paramaters similar to MATLAB
function cgSlow{Tval}(mat, b::Array{Tval,1};
        tol::Real=1e-6, maxits=Inf, maxtime=Inf, verbose=false)

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
                println("CG stopped at maxtime.")
            end
            break
        end
      end

    if verbose
        println("CG stopped after: ", itcnt, " iterations with relative error ", (norm(r)/norm(b)), ".")
    end

    return x
end

#====================================
   PCG Code
=====================================#


# uses BLAS.  As fast as Matlab's pcg.
# set to use similar paramaters
function pcgBLAS{Tval}(mat, b::Array{Tval,1}, pre;
        tol::Real=1e-6, maxits=Inf, maxtime=Inf, verbose=false)

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



        # p = z + beta*p
        BLAS.scal!(n,beta,p,1) # p *= beta
        BLAS.axpy!(1.0,z,p) # p += z

        if (time() - t1) > maxtime
            if verbose
                println("PCG stopped at maxtime.")
            end
            break
        end

       
      end

    if verbose
        println("PCG stopped after: ", itcnt, " iterations with relative error ", (norm(r)/norm(b)), ".")
    end

    return x
end

function pcgSlow{Tval}(mat, b::Array{Tval,1}, pre;
        tol::Real=1e-6, maxits=Inf, maxtime=Inf, verbose=false)


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
                println("PCG stopped at maxtime.")
            end
            break
        end
       
    end

    if verbose
        println("PCG stopped after: ", itcnt, " iterations with relative error ", (norm(r)/norm(b)), ".")
    end

    return x
end
