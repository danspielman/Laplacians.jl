#=
	A data structure that should store triplets of type (neighbor,weight,multiedgeCount) faster
	than an array of arrays of tuples. The backbone of this datastructure will be n linked lists.
	Three methods are implemented:
		- init
		- add
		- purge
=#

## Initially, this will be an array of arrays

type LinkedListStorage{Tv,Ti}
	val::Array{Array{Tuple{Tv,Ti,Ti},1},1}
	n::Ti
end


function llsInit{Tv,Ti}(a::SparseMatrixCSC{Tv,Ti})

	n = a.n

	v = Array{Tuple{Tv,Ti,Ti},1}[]
	for i in 1:n
		push!(v, Tuple{Tv,Ti,Ti}[])
	end

	return LinkedListStorage(v, n)

end


function llsAdd{Tv,Ti}(lls::LinkedListStorage{Tv,Ti}, pos::Ti, t::Tuple{Tv,Ti,Ti})
	push!(lls.val[pos], t)
end


function llsPurge{Tv,Ti}(lls::LinkedListStorage{Tv,Ti}, pos::Ti, rho::Ti, auxVal::Array{Tv,1}, auxMult::Array{Ti,1})

	col = pos
	v = lls.val[pos]

	multSum::Ti = 0
    diag::Tv = 0
    for i in 1:length(v)
        neigh = v[i][3]
        w = v[i][1]
        e = v[i][2]

        auxVal[neigh] += w
        diag += w
        auxMult[neigh] += e
    end

    res = Tv[]
    mult = Ti[]
    ind = Ti[]
    
    for i in 1:length(v)
        neigh = v[i][3]

        if auxVal[neigh] != 0
            @assert(col < neigh, "col < neigh in purge")
            if neigh == col
                @assert(false, "col = neigh in purge")
            else
                # we want to cap the number of multiedges at rho
                actualMult = min(rho, auxMult[neigh])

                push!(res, auxVal[neigh])
                push!(mult, actualMult)
                push!(ind, neigh)

                multSum = multSum + actualMult

                auxVal[neigh] = 0
                auxMult[neigh] = 0
            end
        end
    end

    return diag, res, multSum, mult, ind

end