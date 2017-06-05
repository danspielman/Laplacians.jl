#=

Code for transforming between Laplacian and SDDM solvers.

By Dan Spielman

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
    t1 = time()
    sol = solver(a)
    if verbose
        println("Solver build time: ", round((time() - t1),3), " seconds.")
    end

    f = function(b; verbose=false, jnkargs...)
        if length(pcgIts) > 0
            pcgIts[1] = 1
        end

        t1 = time()
        if isa(sol,Factorization)
            x = sol \ b
        else
            x = sol(b; params...)
        end

        if verbose
            println("Solve time: ", round((time() - t1),3), " seconds.")
        end

        return x
    end
    return f
end

function wrapInterface(solver::Function)
    f = function(a::AbstractArray; tol=0, maxits=Inf, maxtime=Inf, verbose=false, pcgIts=Int[], params...)
        return wrapInterface(solver, a; pcgIts=Int[], verbose=verbose, params...)
    end
    return f
end

function nullSolver(a;params...)
    return 0.0
end

"""
    solveSDDM = cholSDDM(sddm::AbstractMatrix; tol, maxits, maxtime, verbose, pcgIts=Int[])

This functions wraps cholfact so that it satsfies our interface.
It ignores all the keyword arguments.
"""
cholSDDM = wrapInterface(cholfact)


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
    f = lapWrapConnected(sddmSolver, a::AbstractMatrix; kwargs...)

Applies a `sddmSolver` to the Laplacian of the adjacency matrix `a` of a connected graph.
Passes on kwargs to the solver.
`sddmSolver` should be a solver that obeys the interface.
"""
function lapWrapConnected(solver, a::AbstractMatrix; tol::Real=1e-6, maxits=Inf, maxtime=Inf, verbose=false, pcgIts=Int[], params...)
    la = forceLap(a)
    N = size(la)[1]

    ind = findmax(diag(la))[2]
    leave = [1:(ind-1);(ind+1):N]

    lasub = la[leave,leave]
    subSolver = solver(lasub; tol=tol, maxits=maxits, maxtime=maxtime, verbose=verbose, pcgIts=pcgIts);

    tol_=tol
    maxits_=maxits
    maxtime_=maxtime
    verbose_=verbose
    pcgIts_=pcgIts

    f = function(b; tol=tol_, maxits=maxits_, maxtime=maxtime_, verbose=verbose_, pcgIts=pcgIts_)

        bs = b[leave] - mean(b)

        xs = subSolver(bs, tol=tol, maxits=maxits, maxtime=maxtime, verbose=verbose, pcgIts=pcgIts)

        x = zeros(b)
        x[leave] = xs
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

        if length(pcgIts) > 0
            pcgIts[1] = 0
            pcgTmp = Int[0]
        else
            pcgTmp = Int[]
        end


        x = zeros(size(b))
        for i in 1:length(comps)
            ind = comps[i]
            bi = b[ind]
            x[ind] = (solvers[i])(bi;  tol=tol, maxits=maxits, maxtime=maxtime, verbose=verbose, pcgIts=pcgTmp)
            if length(pcgIts) > 0
                pcgIts[1] = max(pcgIts[1],pcgTmp[1])
            end

        end
        return x
    end

end



"""
    f = lapWrapComponents(solver, a::AbstractArray; tol::Real=1e-6, maxits=Inf, maxtime=Inf, verbose=false, pcgIts=Int[], params...)

Applies a Laplacian `solver` that satisfies our interface to each connected component of the graph with adjacency matrix `a`.
Passes kwargs on the solver.
"""
function lapWrapComponents(solver, a::AbstractArray; tol::Real=1e-6, maxits=Inf, maxtime=Inf, verbose=false, pcgIts=Int[], params...)

    t1 = time()

    co = components(a)

    if maximum(co) == 1

        s = solver(a; tol=tol, maxits=maxits, maxtime=maxtime, verbose=verbose, pcgIts=pcgIts, params... )
        if verbose
            println("Solver build time: ", round((time() - t1),3), " seconds.")
        end

        # f(b; tol=tol_, maxits=maxits_, maxtime=maxtime_, verbose=verbose_, pcgIts=pcgIts_) =
        return s

    else

        comps = vecToComps(co)

        solvers = []
        for i in 1:length(comps)
            ind = comps[i]

            asub = a[ind,ind]

            if (length(ind) == 1)
                subSolver = nullSolver

            elseif (length(ind) < 50)
                subSolver = lapWrapConnected(cholSDDM,asub)

            else

                subSolver = solver(asub; tol=tol, maxits=maxits, maxtime=maxtime, verbose=verbose, pcgIts=pcgIts, params... );

            end
            push!(solvers, subSolver)
        end

        if verbose
            println("Solver build time: ", round((time() - t1),3), " seconds.")
        end

        return blockSolver(comps,solvers; tol=tol, maxits=maxits, maxtime=maxtime, verbose=verbose, pcgIts=pcgIts)

    end
end


function lapWrapComponents(solver::Function)
    f(a::AbstractArray; tol=1e-6, maxits=Inf, maxtime=Inf, verbose=false, pcgIts=Int[], params...) = lapWrapComponents(solver, a; tol=tol, maxits=maxits, maxtime=maxtime, verbose=verbose, pcgIts=pcgIts, params... )
    return f
end



"""
    f = lapWrapSDDM(sddmSolver, A::AbstractArray; tol::Real=1e-6, maxits=Inf, maxtime=Inf, verbose=false, pcgIts=Int[], params...)
    f = lapWrapSDDM(sddmSolver)

Uses a `sddmSolver` to solve systems of linear equations in Laplacian matrices.
"""
function lapWrapSDDM(sddmSolver, a::AbstractArray; tol::Real=1e-6, maxits=Inf, maxtime=Inf, verbose=false, pcgIts=Int[], params...)

    @assert(minimum(a) >= 0, "A must be nonnegative")
    @assert(sum(abs(diag(a))) == 0, "A must have zero diagonal")

    f = Laplacians.lapWrapComponents(Laplacians.lapWrapConnected(sddmSolver))(a, tol=tol, maxits=maxits, maxtime=maxtime, verbose=verbose, pcgIts=pcgIts)

    return f

end

function lapWrapSDDM(sddmSolver)

    return Laplacians.lapWrapComponents(Laplacians.lapWrapConnected(sddmSolver))

end



"""
    solver = cholLap(A::AbstractArray)

Uses Cholesky Factorization to solve systems in Laplacians.
"""
cholLap = lapWrapSDDM(cholSDDM)



"""
    f = sddmWrapLap(lapSolver, sddm::AbstractArray; tol::Real=1e-6, maxits=Inf, maxtime=Inf, verbose=false, pcgIts=Int[], params...)

Uses a `lapSolver` to solve systems of linear equations in sddm matrices.
"""
function sddmWrapLap(lapSolver, sddm::AbstractArray; tol::Real=1e-6, maxits=Inf, maxtime=Inf, verbose=false, pcgIts=Int[], params...)

    # Make a new adj matrix, a1, with an extra entry at the end.
    a, d = adj(sddm)
    a1 = extendMatrix(a,d)
    F = lapSolver(a1, tol=tol, maxits=maxits, maxtime=maxtime, verbose=verbose, pcgIts=pcgIts, params...)

    # make a function that solves the extended system, modulo the last entry
    tol_=tol
    maxits_=maxits
    maxtime_=maxtime
    verbose_=verbose
    pcgIts_=pcgIts

    f = function(b; tol=tol_, maxits=maxits_, maxtime=maxtime_, verbose=verbose_, pcgIts=pcgIts_)
        @time xaug = F([b; -sum(b)], tol=tol, maxits=maxits, maxtime=maxtime, verbose=verbose, pcgIts=pcgIts)
        xaug = xaug - xaug[end]
        return xaug[1:a.n]
    end

    return f

end

"""
    f = wrapCaptureRhs(sola::Function, rhss; tol::Real=1e-6, maxits=Inf, maxtime=Inf, verbose=false, pcgIts=Int[], params...)

Captures all the right-hand-sides that are passed to the solver `sola`.  It pushes them into an array called rhhs.
For example

```julia
julia> rhss = []
julia> a = wtedChimera(100)
julia> sola = approxCholLap(a)
julia> wrappedSolver = wrapCaptureRhs(sola,rhss)
julia> b = randn(100)
julia> x = wrappedSolver(b,verbose=true)

PCG BLAS stopped after: 0.0 seconds and 11 iterations with relative error 3.160275810360986e-7.

julia> length(rhss[1])

100
```

"""
function wrapCaptureRhs(sola::Function, rhss; tol::Real=1e-6, maxits=Inf, maxtime=Inf, verbose=false, pcgIts=Int[], params...)

    tol_=tol
    maxits_=maxits
    maxtime_=maxtime
    verbose_=verbose
    pcgIts_=pcgIts

    f = function(b::AbstractArray; tol=tol_, maxits=maxits_, maxtime=maxtime_, verbose=verbose_, pcgIts=pcgIts_)
        push!(rhss, b)
        x = sola(b, tol=tol, maxits=maxits, maxtime=maxtime, verbose=verbose, pcgIts=pcgIts)
    end
end

"""
    f = wrapCapture(solver::Function, mats, rhss; tol::Real=1e-6, maxits=Inf, maxtime=Inf, verbose=false, pcgIts=Int[], params...)

This wraps a solver so that we can capture all the matrices that it solves and all the right-hand-sides.
Those are pushed into the arrays `mats` and `rhss`.
For example

```julia
julia> mats = []
julia> rhss = []
julia> solver = wrapCapture(approxCholLap, mats, rhss)
julia> a = chimera(10)
julia> f = solver(a);
julia> size(mats[1])
(10,10)
julia> b = randn(10)
julia> x = f(b);
julia> rhss
1-element Array{Any,1}:
 [0.404962,-0.827718,0.704616,-0.403223,0.204891,-0.505589,0.907015,1.90266,-0.438115,0.0464351]
```
"""

function wrapCapture(solver::Function, mats, rhss)
    f = function(a::AbstractArray; tol::Real=1e-6, maxits=Inf, maxtime=Inf, verbose=false, pcgIts=Int[], params...)
        tol_=tol
        maxits_=maxits
        maxtime_=maxtime
        verbose_=verbose
        pcgIts_=pcgIts

        push!(mats,a)
        sol1 = solver(a, tol=tol_, maxits=maxits_, maxtime=maxtime_, verbose=verbose_, pcgIts=pcgIts_, params...)
        return wrapCaptureRhs(sol1, rhss, tol=tol_, maxits=maxits_, maxtime=maxtime_, verbose=verbose_, pcgIts=pcgIts_)
    end
    return f
end
