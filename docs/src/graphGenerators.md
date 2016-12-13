# Generators

`Laplacians.jl` implements generators for many standard graphs.
The [`chimera`](@ref) and [`wtedChimera`](@ref) generators are
designed to stress code by combining these standard graphs in tricky
ways.  While no one of these graphs need be a hard case for any
application, the goal is for these generators to explore the space of
graphs in such a way that running on many of them should exercise your
code.

`chimera(n)` generates a random chimera graph.
`chimera(n,k)` first sets the seed of the psrg to k.
In this way, it generates the kth chimera graph, and messes with your
psrg.
`wtedChimera` is similar, but it generates weighted graphs.

## Function list

```@index
Order = [:type, :function]
Pages   = ["graphGenerators.md"]
```

```@autodocs
Modules = [Laplacians]
Pages   = ["graphGenerators.jl"]
Private = false
```
