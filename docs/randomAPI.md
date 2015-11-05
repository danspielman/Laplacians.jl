# random
### randperm
```rst
..  randperm([rng,] n)

Construct a random permutation of length ``n``. The optional ``rng`` argument
specifies a random number generator, see :ref:`Random Numbers <random-numbers>`.
```

Randomly permutes the vertex indices


```julia
randperm(r::AbstractRNG, n::Integer)
randperm(n::Integer)
randperm(mat::AbstractArray{T,2})
randperm(f::Expr)
```

randperm(r::AbstractRNG, n::Integer) at random.jl:1341



