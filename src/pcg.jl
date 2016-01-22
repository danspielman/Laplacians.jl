#=

Implementations of cg and pcg that use BLAS when they can.

Started by Dan Spielman.
Contributors:

=#

# KMP_LOGGING=false
# KMP_PCGITERS=0

# export KMP_LOGGING
# export KMP_PCGITERS


#=
"""
`cg(mat, b; tol, maxits)` solves a symmetric linear system.
`tol` is set to 1e-6 by default,
`maxits` defaults to 100
"""
function cg end

"""
`pcg(mat, b, pre; tol, maxits)` solves a symmetric linear system using preconditioner `pre`.
`pre` should be a function
`tol` is set to 1e-6 by default,
`maxits` defaults to 100
"""
function pcg end
=#

function cg(mat, b::Array{Float64,1}; tol::Real=1e-6, maxits::Integer=100)
    cgBLAS(mat, b, tol=tol, maxits=maxits)
end

function cg(mat, b::Array{Float32,1}; tol::Real=1e-6, maxits::Integer=100)
    cgBLAS(mat, b, tol=tol, maxits=maxits)
end

function cg(mat, b; tol::Real=1e-6, maxits::Integer=100)
    cgSlow(mat, b, tol=tol, maxits=maxits)
end


function pcg(mat, b::Array{Float64,1}, pre; tol::Real=1e-6, maxits::Integer=100)
    pcgBLAS(mat, b, pre, tol=tol, maxits=maxits)
end

function pcg(mat, b::Array{Float32,1}, pre; tol::Real=1e-6, maxits::Integer=100)
    pcgBLAS(mat, b, pre, tol=tol, maxits=maxits)
end

function pcg(mat, b, pre; tol::Real=1e-6, maxits::Integer=100)
    pcgSlow(mat, b, pre, tol=tol, maxits=maxits)
end


# uses BLAS.  As fast as Matlab's pcg.
# set to use similar paramaters
function cgBLAS{Tval}(mat, b::Array{Tval,1};
        tol::Real=1e-6, maxits::Integer=100)

    n = size(mat,2)
    
    x = zeros(Tval,n)
    
    tol = tol * norm(b)
    r = copy(b) 
    p = copy(r)

    # the following would have higher accuracy
    #       rho = sum(r.^2)
    rho = norm(r)^2
    
    for iter in 1:maxits
    
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

       
      end
    return x
end

# For when you can't use BLAS.
# set to use similar paramaters similar to MATLAB
function cgSlow{Tval}(mat, b::Array{Tval,1};
        tol::Real=1e-6, maxits::Integer=100)

    n = size(mat,2)
    
    x = zeros(Tval,n)
    
    tol = tol * norm(b)
    r = copy(b) 
    p = copy(r)

    # the following would have higher accuracy
    #       rho = sum(r.^2)
    rho = norm(r)^2
    
    for iter in 1:maxits
    
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
       
      end

    return x
end

#====================================
   PCG Code
=====================================#


# uses BLAS.  As fast as Matlab's pcg.
# set to use similar paramaters
function pcgBLAS{Tval}(mat, b::Array{Tval,1}, pre;
        tol::Real=1e-6, maxits::Integer=100)

    n = size(mat,2)
    
    x = zeros(Tval,n)
    
    tol = tol * norm(b)
    r = copy(b) 
    z = pre(r)
    p = copy(z)

    rho = dot(r, z) # BLAS.dot does not seem faster

    itcnt = 0
    for iter in 1:maxits
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

       
      end

    #=
    if KMP_LOGGING
        KMP_PCGITERS = itcnt
        println("ITERS : ", itcnt)
    end

    println(KMP_LOGGING)
    println("iters : ", itcnt)
    =#
    
    return x
end

function pcgSlow{Tval}(mat, b::Array{Tval,1}, pre;
        tol::Real=1e-6, maxits::Integer=100)


    n = size(mat,2)
    
    x = zeros(Tval,n)
    
    tol = tol * norm(b)
    r = copy(b) 
    z = pre(r)
    p = copy(z)

    rho = dot(r, z) 
    
    for iter in 1:maxits
    
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
       
      end
    return x
end
