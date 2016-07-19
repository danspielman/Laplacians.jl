#=
	A function that returns an epsilon sparsifier of a graph.
=#

include("condNumber.jl")

function sparsify{Tv,Ti}(a::SparseMatrixCSC{Tv,Ti}, eps::Tv; treeAlg=akpw, cnTol=1e-5, verbose=false)
	
	n = a.n;

	tree = treeAlg(a)
	st = compStretches(tree, a - tree)

	(ai,aj,av) = findnz(triu(a - tree))
	(si,sj,sv) = findnz(triu(st))

	#=
		Next, we will guess the number of edges to add back to the tree. 
		We know from Spielman Strivsastava that we can get an eps sparsifier with eps^-2 n log n edges.
		So, that quantity should be a good upper bound.
		Given that most of the graphs in chimera are already sparse, we can try generating an epsilon
		sparsifier with fewer edges than that. 
	=#

	# k = ceil(Int64, eps^-2 * n * log(n));
	k = min(ceil(Int64, qn), ceil(Int64, log(nnz(a))))

	println(k, " ", nnz(a))

	# if the graph is sparse enough, an epsilon aproximation with few edges can be the initial graph
	if k >= nnz(a)
		return a
	end

	while true

		if verbose
			println("Trying a sparsifier with ", k + (n - 1), " edges")
		end

		probs = sv / sum(sv)
		take = (k * rand(length(probs))) .> probs

		sp = sparse(ai[take], aj[take], av[take], a.n, a.m)
		sp = sp + sp' + tree

		cn = condNumber(lap(a), lap(sp), tol=cnTol)

		if verbose
			println("Got condition number between ", cn)
		end

		if cn[2] < (1 + eps) / (1 - eps)
			return sp
		else
			k = ceil(Int64, k * 1.1)
		end

	end

end