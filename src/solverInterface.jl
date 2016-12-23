#=

Code for solving Laplacian and Diagnally Dominant Systems
This is not all of it, just the short and sweet routines.

Started by Dan Spielman

  lapWrapSolver: takes a solver for DD systems, and uses it to solve a lap system in la
  lapWrapSolver(solver, la::AbstractArray)
  lapWrapSolver(solver)
  lapWrapSolver(solver, la::AbstractArray, b) = lapWrapSolver(solver,la)(b)    

  For example, to make a Cholesky-based solver for Laplacians, we created
  lapChol = lapWrapSolver(cholfact)


  augmentTree : takes a spanning tree, a graph, and adds in 2k edges of high stretch
    takes both as sparse matrices


  augTreeSolver : a solver that shouldn't suck

=#


"""
    solveA = wrapInterface(solver::Function, A::AbstractMatrix; tol, maxits, maxtime, verbose, pcgIts=Int[],params...)
    solverConstructor = wrapInterface(A::AbstractMatrix; tol, maxits, maxtime, verbose, pcgIts=Int[],params...)

Returns a function that discards `tol`, `maxits`, `maxtime` and `verbose`, 
sets `pcgIts` to 0 (because it might not be using pcg), 
and passes whatever `params` are left to the solver.

# Examples

```julia
julia> a = randn(5,5);
julia> a = a * a';
julia> solvea = wrapInterface(cholfact, a, maxits=100, verbose=true);
julia> b = randn(5,1);
julia> norm(a*solvea(b, verbose=false)-b)
1.575705319704736e-14

julia> f = wrapInterface(cholfact)
julia> solvea = f(a, maxits=1000, maxtime = 1)
julia> norm(a*solvea(b, verbose=false, maxtime = 10)-b)
1.575705319704736e-14
```
"""
function wrapInterface(solver::Function, a::AbstractMatrix; tol=0, maxits=Inf, maxtime=Inf, verbose=false, pcgIts=Int[],params...)
    sol = solver(a)
    f = function(b; jnk=false, jnkargs...)
        if length(pcgIts) > 0
            pcgIts[1] = 1
        end
        if isa(sol,Factorization)
            x = sol \ b
        else
            x = sol(b; params...)
        end
        return x
    end
    return f
end

function wrapInterface(solver::Function)
    f = function(a::AbstractArray; tol=0, maxits=Inf, maxtime=Inf, verbose=false, pcgIts=Int[], params...)
        return wrapInterface(solver, a; pcgIts=Int[], params...)
    end
    return f
end


"""
    solveSDD = cholSDD(SDD::AbstractMatrix; tol, maxits, maxtime, verbose, pcgIts=Int[])

This functions wraps cholfact so that it satsfies our interface.
It ignores all the keyword arguments.
"""    
cholSDD = Laplacians.wrapInterface(cholfact)


"""
    la = forceLap(a)

Create a Laplacian matrix from an adjacency matrix. 
If the input looks like a Laplacian, throw a warning and convert it.
"""
function forceLap(a::AbstractArray)

    if minimum(a) < 0
        warn("The input should be an adjacency matrix, whereas this one has negative entries.")
        af = abs(a)
        af = af - spdiagm(diag(af))
    elseif sum(abs(diag(a))) > 0
        warn("The input should be an adjacency matrix, whereas this one has diagonal entries.")
        af = a - spdiagm(diag(a))
    else
        af = a
    end

    return spdiagm(vec(sum(af,1))) - af
end


"""
    f = lapWrapConnected(sddSolver, a::AbstractMatrix; kwargs...)

Applies a `sddSolver` to the Laplacian of the adjacency matrix `a` of a connected graph.
Passes on kwargs to the solver.
`sddSolver` should be a solver that obeys the interface.
"""
function lapWrapConnected(solver, a::AbstractMatrix; tol::Real=1e-6, maxits=Inf, maxtime=Inf, verbose=false, pcgIts=Int[], params...)
    la = forceLap(a)
    N = size(la)[1]
    lasub = la[1:(N-1),1:(N-1)]
    subSolver = solver(lasub; tol=tol, maxits=maxits, maxtime=maxtime, verbose=verbose, pcgIts=pcgIts);

    tol_=tol
    maxits_=maxits
    maxtime_=maxtime
    verbose_=verbose
    pcgIts_=pcgIts

    f = function(b; tol=tol_, maxits=maxits_, maxtime=maxtime_, verbose=verbose_, pcgIts=pcgIts_)

        bs = b[1:(N-1)]
        
        xs = subSolver(bs, tol=tol, maxits=maxits, maxtime=maxtime, verbose=verbose, pcgIts=pcgIts)

        x = [xs;0]
        x = x - mean(x)
        return x
    end
    
    return f
end


function lapWrapConnected(solver::Function)
    f(a::AbstractArray; tol=1e-6, maxits=Inf, maxtime=Inf, verbose=false, pcgIts=Int[], params...) = lapWrapConnected(solver, a; tol=tol, maxits=maxits, maxtime=maxtime, verbose=verbose, pcgIts=pcgIts, params... )
    return f
end
        

"""Apply the ith solver on the ith component"""
function blockSolver(comps, solvers; tol::Real=1e-6, maxits=Inf, maxtime=Inf, verbose=false, pcgIts=Int[])

    tol_=tol
    maxits_=maxits
    maxtime_=maxtime
    verbose_=verbose
    pcgIts_=pcgIts

    f = function(b; tol=tol_, maxits=maxits_, maxtime=maxtime_, verbose=verbose_, pcgIts=pcgIts_)
        x = zeros(size(b))
        for i in 1:length(comps)
            ind = comps[i]
            bi = b[ind]
            x[ind] = (solvers[i])(bi;  tol=tol, maxits=maxits, maxtime=maxtime, verbose=verbose, pcgIts=pcgIts)
        end
        return x
    end
        
end



"""
    f = lapWrapComponents(solver, a::AbstractArray; kwargs)

Applies a Laplacian `solver` that satisfies our interface to each connected component of the graph with adjacency matrix `a`.
Passes kwargs on the solver.
"""
function lapWrapComponents(solver, a::AbstractArray; tol::Real=1e-6, maxits=Inf, maxtime=Inf, verbose=false, pcgIts=Int[], params...)

    co = components(a)

    if maximum(co) == 1
        # f(b; tol=tol_, maxits=maxits_, maxtime=maxtime_, verbose=verbose_, pcgIts=pcgIts_) =         
        return solver(a; tol=tol, maxits=maxits, maxtime=maxtime, verbose=verbose, pcgIts=pcgIts, params... )

    else
        
        comps = vecToComps(co)

        solvers = []
        for i in 1:length(comps)
            ind = comps[i]
            
            asub = a[ind,ind]

            if (length(ind) == 1)
                ssubSolver = x->0
            
            elseif (length(ind) < 50)
                subSolver = lapWrapConnected(cholfact,asub)

            else

                subSolver = solver(asub; tol=tol, maxits=maxits, maxtime=maxtime, verbose=verbose, pcgIts=pcgIts, params... );

            end
            push!(solvers, subSolver)
        end

        return blockSolver(comps,solvers; tol=tol, maxits=maxits, maxtime=maxtime, verbose=verbose, pcgIts=pcgIts)

    end
end


function lapWrapComponents(solver::Function)
    f(a::AbstractArray; tol=1e-6, maxits=Inf, maxtime=Inf, verbose=false, pcgIts=Int[], params...) = lapWrapComponents(solver, a; tol=tol, maxits=maxits, maxtime=maxtime, verbose=verbose, pcgIts=pcgIts, params... )
    return f
end

    

"""
    f = lapWrapSDD(sddmSolver, A::AbstractArray; tol::Real=1e-6, maxits=Inf, maxtime=Inf, verbose=false, pcgIts=Int[], params...)
    f = lapWrapSDD(sddmSolver)

Uses a `sddmSolver` to solve systems of linear equations in Laplacian matrices.
"""
function lapWrapSDD(sddmSolver, a::AbstractArray; tol::Real=1e-6, maxits=Inf, maxtime=Inf, verbose=false, pcgIts=Int[], params...)

    @assert(minimum(a) >= 0, "A must be nonnegative")
    @assert(sum(abs(diag(a))) == 0, "A must have zero diagonal")

    f = Laplacians.lapWrapComponents(Laplacians.lapWrapConnected(sddmSolver))(a, tol=tol, maxits=maxits, maxtime=maxtime, verbose=verbose, pcgIts=pcgIts)

    return f
                                     
end

function lapWrapSDD(sddmSolver)

    return Laplacians.lapWrapComponents(Laplacians.lapWrapConnected(sddmSolver))

end



"""
    solver = cholLap(A::AbstractArray)

Uses Cholesky Factorization to solve systems in Laplacians.
"""    
cholLap = lapWrapSDD(cholSDD)
    


