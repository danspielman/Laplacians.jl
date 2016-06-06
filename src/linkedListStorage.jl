#=
	A data structure that should store triplets of type (neighbor,weight,multiedgeCount) faster
	than an array of arrays of tuples. The backbone of this datastructure will be n linked lists.
	Three methods are implemented:
		- init
		- add
		- purge

	Each element stores the following information:
		- edge weight
		- multiedge count
		- neighboring vertex
		- position of next element in LinkedListStorage
=#

type element{Tv,Ti}
	edgeWeight::Tv
	edgeCount::Ti
	neighbor::Ti
	next::Ti
end

type LinkedListStorage{Tv,Ti}
	val::Array{element{Tv,Ti},1}		# a big block of memory storing the values from the linked lists

	free::Array{Ti,1}					# an array storing all the positions of free memory in val. it is circular
	left::Ti
	right::Ti

	first::Array{Ti,1}					# the position of the first element in i's linked list
	last::Array{Ti,1}					# the position of the last element in i's linked list

	n::Ti 								# the number of vertices in the graph
	size::Ti 							# total size of the data structure
end

function llsInit{Tv,Ti}(a::SparseMatrixCSC{Tv,Ti}, rho::Ti)

	n = a.n
	m = length(a.nzval)
	size = rho * m

	first = -ones(Ti, n)
	last = -ones(Ti, n)

	val = Array{element{Tv,Ti},1}(size)

	free = collect(Ti, 1:size)
	left = right = size

	return LinkedListStorage(val, free, left, right, first, last, n, size)

end

function llsAdd{Tv,Ti}(lls::LinkedListStorage{Tv,Ti}, pos::Ti, t::Tuple{Tv,Ti,Ti})

	if lls.last[pos] == -1
		lls.left = moduloNext(lls.left, lls.size)
		ind = lls.free[lls.left];

		lls.first[pos] = lls.last[pos] = ind
		lls.val[ind] = element(t[1], t[2], t[3], -1)
	else
		lls.left = moduloNext(lls.left, lls.size)
		ind = lls.free[lls.left];

		# set the new value for next
		lls.val[lls.last[pos]].next = ind

		lls.last[pos] = ind
		lls.val[ind] = element(t[1], t[2], t[3], -1)
	end

end


function llsPurge{Tv,Ti}(lls::LinkedListStorage{Tv,Ti}, pos::Ti, auxVal::Array{Tv,1}, auxMult::Array{Ti,1})

 	multSum::Ti = 0
 	diag::Tv = 0

	i = lls.first[pos]
	while i != -1
		neigh = lls.val[i].neighbor
		w = lls.val[i].edgeWeight
		e = lls.val[i].edgeCount

		auxVal[neigh] += w
		diag += w
		auxMult[neigh] += e

		i = lls.val[i].next
	end

    res = Tv[]
    mult = Ti[]
    ind = Ti[]

    i = lls.first[pos]
    while i != -1
        neigh = lls.val[i].neighbor

        if auxVal[neigh] != 0
            @assert(pos < neigh, "current element < neigh in purge")
            if neigh == pos
                @assert(false, "current element = neigh in purge")
            else
                # TODO: we want to cap the number of multiedges
                actualMult = auxMult[neigh]

                push!(res, auxVal[neigh])
                push!(mult, actualMult)
                push!(ind, neigh)

                multSum = multSum + actualMult

                auxVal[neigh] = 0
                auxMult[neigh] = 0
            end
        end

        # i was just freed
        lls.right = moduloNext(lls.right, lls.size)
        lls.free[lls.right] = i

        i = lls.val[i].next
    end

    return diag, res, multSum, mult, ind

end

function moduloNext{Ti}(i::Ti, maxSize::Ti)
	return i < maxSize ? i + 1 : 1
end