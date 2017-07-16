#=
Written by Rasmus Kyng.
Based closely on Matlab code by Anup Rao, Sushant Sachdeva & Rasmus Kyng.
Matlab code available at https://github.com/sachdevasushant/Isotonic

Interior point method for computing the Isotonic regression (under Euclidean norm) on arbitrary directed acyclic graphs. The code accompanies the paper 'Fast, Provable Algorithms for Isotonic Regression in all l_p-norms' at NIPS 2015

This repository presently contains code for computing isotonic regression in directed acyclic graphs (DAG).

isotonicIPM takes as input the adjacency matrix ‘A’ of a DAG, followed by a vector ‘v’ giving the initial values on the vertices. v should be a vector of length equal to the number of vertices.

On input ‘A’ and ‘v’, the code computes a vector ‘x’ of same length as v, such that x minimizes the quantity: sum_i (x(i) - v(i))^2 subject to x(i) <= x(j) if (i,j) is an entry in ‘a’.

=#


"""Compute isotonic regression of v with constraints given by A,
and multiplicative error param eps,
using a PD SDD linear system solver, given by argument *solver*."""
function isotonicIPMrelEps{Tv,Ti}(A::SparseMatrixCSC{Tv,Ti},
                                  v::Array{Tv,1},
                                  eps::Real=0.1,
                                  solver=(H -> augTreeSddm(H,tol=1e-1,maxits=1000)))
    return isotonicIPM(A,v,eps,solver,true)
end

"""Compute isotonic regression of v with constraints given by A,
and additive error param eps,
using a PD SDD linear system solver, given by argument *solver*.
Setting relGapTermination=true means the IPM continues until a multiplicative error of eps is reached."""
function isotonicIPM{Tv,Ti}(A::SparseMatrixCSC{Tv,Ti},
                            v::Array{Tv,1},
                            eps::Tv=0.1,
                            solver=(H -> augTreeSddm(H,tol=1e-1,maxits=1000)),
                            relGapTermination::Bool=false)
    n = size(A)[1]
    topoOrder = toposort(A)
    permMat = (speye(n))[1:n,topoOrder]
    Bt = -dirEdgeVertexMat(A) #sign convention implies all edges should get positive values
    m = size(Bt)[1]
    
    Theta = m
    eps0 = eps
    beta2 = 10.0

    x = 0.0+collect(1:n) # our initial guess
    x = permMat*x #this seems to be the right direction for permMat
    
    if notFeasible(x, Bt)
        error("Input graph is not a DAG!")
    end

    eTaStartFn(y) = ( (y-v)'*Bt'*(1./(Bt*y))/norm(y-v)^2 )[1]
    #This seems to be the formula for the minimizer?

    #TODO
    #eTaStartFn(y) = ( (y-v)'*Bt'*(1./(Bt*y))/norm(y-v) )[1]
    #This is the formula Anup used?
    #Seems to give better results?
    
    eTa1 = eTaStartFn(x)
    #@printf "eTa1 = %f\n" eTa1
    while eTa1<0
        x = 2*x
        eTa1 = eTaStartFn(x)
        #Note: y'*Bt'*(1./(Bt*y)) is always pos when y is isotonic,
        #We can show this loop will complete
        #println("doubling")
    end

    #@printf "eTa1 = %f\n" eTa1
        
    if eTa1 < 0
        error("bad initialization!")
    end
    eTa2 = Theta/eps0
    mu0 = eTa1

    relGap = Inf
    i=0;    
    muTooSmall = mu0 < eTa2
    
    # numSteps = 0

    while muTooSmall
        F(y) = mu0*norm(y-v)^2 - sum(log.(Bt*y))
        centered = false
        while !centered
            # numSteps += 1
            # if numSteps % 30 == 1
            #     println("So far ", numSteps, " steps were made")
            # end


            (H,xNewton) = l2NewtonStep( Bt,x,v,mu0, solver )
            
            s = Bt*x
            gradF = 2*mu0*(x-v) - Bt'*(1./s)
            x = backtrackLineSearch(F,gradF,Bt,xNewton,x)
            centMeasure = (xNewton'*H*xNewton)[1]
            centered = centMeasure < 10.0^-2
            i = i+1
        end
        
        dval = dualVal(Bt,x,v,mu0)
        relGap = (sum((x-v).^2) - dval)/sum((x-v).^2)
        if relGapTermination
            mu0 = mu0*beta2
            if relGap < eps
                break
            end
            #muTooSmall stays false: we only terminate if relGap is good
        else
            #this is the old termination condition:
            #stop if within 1 % relGap or if eta is too large?
            if relGap < 10.0^-2
                break
            end
            muTooSmall = mu0 < eTa2
            if muTooSmall
                mu0 = mu0*beta2
            end
        end
    end
    
    if notFeasible(x, Bt)
        error("output is not isotonic!")
    end
    
    return (x,relGap,i)
end
        
function l2NewtonStep{Tv,Ti}( Bt::SparseMatrixCSC{Tv,Ti},
                              x::Array{Tv,1},
                              v::Array{Tv,1},
                              mu0::Real,
                              solver )
    m = size(Bt)[1]
    n = size(Bt)[2]
    d1 = 2 * mu0 * speye(n)
    s = Bt * x
    d2 = sparse(1:m, 1:m, s .^ -2, m, m)
    grad = 2 * mu0 * (x - v) - Bt'*(1./s)
    H1 = sparse(Bt'*d2*Bt)

    H = H1 + d1

    invMinEntryH = maximum(abs.(s)) ^ 2

    #println("Maximum eigenvalue of the matrix we are computing: ", eigs(invMinEntryH * H;nev=1,which=:LM,tol=1e-1)[1][1])
    #println("Smallest eigenvalues of the matrix we are computing: ", eigs(invMinEntryH * H;nev=2,which=:SM,tol=1e-1)[1][1])
    F = solver(invMinEntryH * H)
    xNewton = -F(grad)*invMinEntryH
    
    return (H,xNewton)
end

function dualVal{Tv,Ti}( Bt::SparseMatrixCSC{Tv,Ti},
                         x::Array{Tv,1},
                         v::Array{Tv,1},
                         mu0::Real)
    s = Bt*x # all pos entries
    f = 1/mu0 * 1./s 
    dval = (-1/4 * sum((Bt'*f).^2) - f'*Bt*v)[1]
    return dval
end

function backtrackLineSearch{Tv,Ti}(F,
                             gradF::Array{Tv,1},
                             Bt::SparseMatrixCSC{Tv,Ti},
                             xNewton::Array{Tv,1},
                             x::Array{Tv,1})
    #track line search as in Boyd's book
    A0 = 0.01;
    B0 = 0.5;

    t = 1;
    iter = 0;

    xt = x + t * xNewton;

    while (notFeasible(xt, Bt)) || (F(xt) > F(x) + A0 * t * (gradF'*xNewton)[1])
        t = B0 * t;
        iter = iter + 1;
        if iter >= 100
            error("Number of iterations in line search exceeded the limit")
        end        
        xt = x + t * xNewton;
    end

    return xt
end

function notFeasible{Tv,Ti}(x::Array{Tv,1},
                     Bt::SparseMatrixCSC{Tv,Ti})
    return minimum(Bt * x) <= 0 #not feasible if some entry is non-positive
end
