# graph Generators

`Laplacians.jl` implements generators for many standard graphs.
It the `chimera` and `wtedChimera` generators combine these standard
graphs by

*  joining disjoint copies with edges,
* forming Kronecker products, and
* taking `generalizedNecklace` products.


```@index
Order = [:type, :function]
Pages   = ["graphGenerators.md"]
```

```@autodocs
Modules = [Laplacians]
Pages   = ["graphGenerators.jl"]
Private = false
```
