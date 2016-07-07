# b is a collection of right hand sides. samplingSolver permits us to solve multiple right hand sides
function pcg(mat, b, pre::Function; tol::Real=1e-6, maxits=Inf, maxtime=Inf, verbose=false)
    pcgSlow(mat, b, pre, tol=tol, maxits=maxits, maxtime=maxtime, verbose=verbose)
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

function matDot{Tv}(a::Array{Tv,2}, b::Array{Tv,2})
end