#=

Primal-dual predictor corrector method  for maxflow
form max_x t s.t. Bf = td, u>=f>=0. 
Input: vertex-edge incidence matrix B,demand d, capacities u 
 
=#


function max_flow_IPM{Tv,Ti}(Bt::SparseMatrixCSC{Tv,Ti},
                            b1::Array{Tv,1},
                            u::Array{Tv,1};
                            lapSolver = (H -> lapWrapSolver(augTreeSolver,H,tol=1e-8,maxits=1000)),
                             tol::Real=1e-6)

    maxIter = 100
    numIter = 0

    m = size(Bt)[2]
    n = size(Bt)[1] 
    A = [Bt spzeros(n,m) -sparsevec(b1); speye(m) speye(m) spzeros(m,1)]
    b = [zeros(n);u]
    c = [zeros(2*m); -1]
    n1 = size(A)[1]
    n2 = size(A)[2]
    trunc = 0.995

    ipm_max_flow_hess_solve = ((x,s) -> ipm_max_flow_shur_solve(Bt,b1,x,s,m,n,lapSolver))

    #computing the intial point
    (x,y,s) = ipm_max_flow_initial_point(A,b,c,u,m,n,ipm_max_flow_hess_solve)
    #display(minimum(s))

    for i = 1:maxIter
     
        numIter = numIter+1

        mu_cur = BLAS.dot(x,s)/n2  #current mu
 
        #predictor step
        rd = (x.*s)
        rb = (A*x - b)
        rc = A'*y + s - c

        norm_prim_feas = norm(rb)
        norm_dual_feas = norm(rc)
        #dual_gap = sum(rd)/m
        @printf("Iteration %d, ||r_p||=%f, ||r_d||=%f, mu=%f\n", numIter, norm_prim_feas, norm_dual_feas, mu_cur);
        if norm_prim_feas <= tol && norm_dual_feas <= tol && mu_cur <= tol
            println("Termination tolerance reached.");
            numIter = i
            break 
        end

        Hinv = ipm_max_flow_hess_solve(x,s)
        (dx,dy,ds) = ipm_max_flow_solver(A,rb,rc,rd,x,s,Hinv)
    
    #compute the step size	
        alpha_p = ipm_max_flow_step_length(x,dx)
        alpha_d = ipm_max_flow_step_length(s,ds)

    #compute corrector step parameters
    	mu_aff = BLAS.dot((x + alpha_p*dx),(s+alpha_d*ds))/n2
    	sig = (min(1,mu_aff/mu_cur))^3
        println("sigma parameter is ",sig)

    #corrector step
    	rd = rd + dx.*ds - (sig*mu_cur)[1,1]
        (dx,dy,ds) = ipm_max_flow_solver(A,rb,rc,rd,x,s,Hinv)
       

    #compute the step size
        alpha_pk = ipm_max_flow_step_length(x,dx)
        alpha_dk = ipm_max_flow_step_length(s,ds)
    

    	x = x + alpha_pk*dx
    	y = y + alpha_dk*dy
    	s = s + alpha_dk*ds

    end

    if numIter >= maxIter
          println("Maximum number of iteration reached.");
    end

    return (x,y,s,numIter)

end

function ipm_max_flow_step_length{Tv}(x::Array{Tv,1},dx::Array{Tv,1},maxstepsize::Float64 = 0.99)

    stepsizes = -x./dx;

    idx_pos = find(stepsizes .> 0);
    if isempty(idx_pos)
        return maxstepsize;
    else
        return minimum([0.999*stepsizes[idx_pos]; maxstepsize]);
    end



end


function ipm_max_flow_initial_point{Tv,Ti}(A::SparseMatrixCSC{Tv,Ti},
                              b::Array{Tv,1},
                              c::Array{Tv,1},
                              u::Array{Tv,1},
                              m,
                              n,
                              ipm_max_flow_hess_solve)


n2 = size(A)[2]
Hinv = ipm_max_flow_hess_solve(ones(n2),ones(n2))
#x0 = A'*Hinv(b)
x0 = -c + A'*Hinv(b + A*c)
#y0 = Hinv(A*c)
y0 = Hinv(A*c - b)
s0 = c - A'*y0


del_x = max(-1.5*minimum(x0),0)
del_s = max(-1.5*minimum(s0),0)
x0 = x0 + del_x
s0 = s0 + del_s

del_x = (0.5*BLAS.dot(x0,s0)/sum(s0))
del_s = 0.5*BLAS.dot(x0,s0)/sum(x0)


x0 = x0 + del_x
s0 = s0 + del_s

#make sure the flow satisfies the capacity constraint
idx_large = find(x0[m+1:m+m] .>= u)
if ~isempty(idx_large)
      x0[n+idx_large] = (2/3)*u[idx_large]
  end

return (x0,y0,s0)

end




function ipm_max_flow_solver{Tv,Ti}(A::SparseMatrixCSC{Tv,Ti},
                        rb::Array{Tv,1},
                        rc::Array{Tv,1},
                        rd::Array{Tv,1},
                        x::Array{Tv,1},
                        s::Array{Tv,1},Hinv)
    

    Sinv = spdiagm((s.^-1)[:,1])
    X = spdiagm(x[:,1])
    rhs = -rb - A*(X*(Sinv*rc)) + A*(Sinv*rd) 
##
    dy = Hinv(rhs)
    ds = -rc - A'*dy
    dx = -Sinv*rd - X*(Sinv*ds) 
##

return(dx,dy,ds)
end



function ipm_max_flow_shur_solve{Tv,Ti}(Bt::SparseMatrixCSC{Tv,Ti},
                          b1::Array{Tv,1},
                          x::Array{Tv,1},
                          s::Array{Tv,1},
                          m::Integer,
                          n::Integer,
                          lapSolver)

    x1 = x[1:m,1]
    x2 = x[m+1:2*m,1]
    x3 = x[2*m+1,1]

    s1 = s[1:m,1]
    s2 = s[m+1:2*m,1]
    s3 = s[2*m+1,1]
    

    d1 = (s1.^-1.*(x1))  #might need to add back the 1e+16 term
    D1 = spdiagm(d1[:,1])

    d2 = (s2.^-1.*(x2))  #might need to add back the 1e+16 term
    D2 = spdiagm(d2[:,1])

    d3 = (s3.^(-1/2))*(x3.^(1/2))

    wt = 1./(1./d1 + 1./d2)
    ratioWt = maximum(wt)/minimum(wt)
    println("Ratio of edge weights: ", ratioWt)

    
    la = Bt*spdiagm(wt[:,1])*Bt' 
   # laInv = lapSolver(la)
    

    dinv = 1./(d1 + d2)
    Dinv = spdiagm(dinv[:,1])


    Sinv = spdiagm((s.^-1)[:,1])
    X = spdiagm(x[:,1])

    function Hinv(rhs)
        rhs1 = rhs[1:n,1]
        rhs2 = rhs[n+1:end,1]

        rr = (rhs1 - Bt*D1*Dinv*rhs2)
        v = d3*b1
     
        #v1 = laInv(v)
        #v2 = laInv(rr)
        u1 = la\v
        u2 = la\rr
        
        df = u2 - BLAS.dot(v,u2)*u1/(1+BLAS.dot(v,u1))
        dw = Dinv*(rhs2 - D1*Bt'*df) 
        dy = [ df;dw]
        return dy 
    end

    return Hinv 

end


