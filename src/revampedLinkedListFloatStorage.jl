
#= 
next will point to the next element
If totalSize % blockSize == 0 for the current position, then next will be pointing to a new block
=#
immutable element{Tv,Ti}
	edgeWeight::Tv
	edgeCount::Tv
	neighbor::Ti
	next::Ti
end

#=
The linked list starts out with a small size, say close to 1000. At each step we increase its size
by sizeIncrease. We hold blockSize consecutive elements in memory for every 
=#
type LinkedListStorage{Tv,Ti}
	val::Array{element{Tv,Ti},1}		# a big block of memory storing the values from the linked lists

	free::Array{Ti,1}					# a circular array storing all the starting positions of blocks of free memory in [left, right)
	left::Ti 							# the first position in the array containing a free element
	right::Ti 							# the end of the sequence of free elements. val[right] is occupied in memory

	first::Array{Ti,1}					# the position of the first element in i's linked list
	last::Array{Ti,1}					# the position of the last element in i's linked list

	size::Ti 							# total size of the data structure (just counts the number of blocks)
	sizeIncrease::Ti                    # the size before each expansion. it should be divisible by blockSize
	blockSize::Ti 						# size of individual blocks of memory
end

function llsInit{Tv,Ti}(a::SparseMatrixCSC{Tv,Ti}; startingSize::Ti = 1000, blockSize::Ti = 20)

	n = a.n
	startingSize::Ti = max(blockSize, startingSize - startingSize % blockSize)
	sizeIncrease::Ti = startingSize

	first = -ones(Ti,n)
	last = -ones(Ti,n)

	val = Array{element{Tv,Ti}}(startingSize * blockSize)

	free = collect(Ti, 1:startingSize)
	left = startingSize
	right = -1

	return LinkedListStorage(val, free, left, right, first, last, startingSize, sizeIncrease, blockSize)

end

function llsAdd{Tv,Ti}(lls::LinkedListStorage{Tv,Ti}, v::Ti, t::Tuple{Tv,Tv,Ti})

	if lls.last[v] == -1
		# nothing was added to the current element's linked list so far

		# if we are out of memory - increase the size of the linked list
		if moduloNext(lls.left, lls.size) == lls.right
			incSize(lls)
		end

		# get the position of a new free block
		lls.left = moduloNext(lls.left, lls.size)
		ind = (lls.free[lls.left] - 1) * lls.blockSize + 1;
		@assert(lls.left != lls.right, "Overwriting values in the linked list struct")

		# this new block will be the starting point for v's linked list
		lls.first[v] = lls.last[v] = ind
		lls.val[ind] = element(t[1], t[2], t[3], -1)

		# set lls.right on the first insertion
		if lls.right == -1
			lls.right = lls.left
		end
	else
		# we're adding element to an already existing linked list

		# check if we still have space in the current block
		if lls.last[v] % lls.blockSize != 0
			ind = lls.last[v] + 1
		else
			# if we are out of memory - increase the size of the linked list
			if moduloNext(lls.left, lls.size) == lls.right
				incSize(lls)
			end

			lls.left = moduloNext(lls.left, lls.size)
			ind = (lls.free[lls.left] - 1) * lls.blockSize + 1;
			@assert(lls.left != lls.right, "Overwriting values in the linked list struct")
		end

		# update the 'next' field for the last element in v' linked list
		prev = lls.last[v]
		lls.val[prev] = element(lls.val[prev].edgeWeight, lls.val[prev].edgeCount, lls.val[prev].neighbor, ind)

		# add in the new element to the linked list
		lls.last[v] = ind
		lls.val[ind] = element(t[1], t[2], t[3], -1)
	end

end

function llsPurge{Tv,Ti}(lls::LinkedListStorage{Tv,Ti}, pos::Ti, auxVal::Array{Tv,1}, auxMult::Array{Tv,1},
	weight::Array{Tv,1}, mult::Array{Tv,1}, ind::Array{Ti,1}; capEdge::Bool=false, rho::Tv=0.0, xhat::Array{Tv,2}=zeros(1,1))

 	multSum::Tv = 0
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

    numPurged = 0

    i = lls.first[pos]
    while i != -1
        neigh = lls.val[i].neighbor

        if auxVal[neigh] != 0
            @assert(pos < neigh, "current element < neigh in purge")
            if neigh == pos
                @assert(false, "current element = neigh in purge")
            else
            	# this is an edge between pos and neigh. if capEdge is true, we will use the effective resistance estimates to cap this edge
                # actualMult = min(auxMult[neigh], rho)
                actualMult = auxMult[neigh]

                if capEdge
                	p = pos
                	q = neigh
                	w = auxVal[neigh]
                	lev = w * norm(xhat[p,:] - xhat[q,:])^2

                	actualMult = min(actualMult, rho * lev)
                end

                numPurged += 1
                weight[numPurged] = auxVal[neigh]
                mult[numPurged] = actualMult
                ind[numPurged] = neigh

                multSum = multSum + actualMult

                auxVal[neigh] = 0
                auxMult[neigh] = 0
            end
        end

        # i was just freed
        if i % lls.blockSize == 0
        	lls.right = moduloNext(lls.right, lls.size)
        	lls.free[lls.right] = i / lls.blockSize
        end

        i = lls.val[i].next
    end

    return diag, multSum, numPurged

end

# increase the size of a lls structure
function incSize{Tv,Ti}(lls::LinkedListStorage{Tv,Ti})

	#=
	we know that if we are increasing the size of the structure, we are out of free blocks
	so, we can add aditional memory at the end of the structure
	=#
	lls.left = lls.size + 1
	lls.right = 1

	append!(lls.val, Array{element{Tv,Ti}}(lls.sizeIncrease * lls.blockSize))
	append!(lls.free, collect(Ti, (lls.size + 1):(lls.size + lls.sizeIncrease)))

	# update the number of free blocks and the number of blocks by which we increase the size every time
	lls.size += lls.sizeIncrease
	lls.sizeIncrease += ceil(Ti, 0.25 * lls.sizeIncrease)
	lls.sizeIncrease = max(lls.blockSize, lls.sizeIncrease - lls.sizeIncrease % lls.blockSize)

end

function moduloNext{Ti}(i::Ti, maxSize::Ti)
	return i < maxSize ? i + 1 : 1
end