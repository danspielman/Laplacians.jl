# IO
### readIJ
to read a simple edge list, each line being an (i, j) pair


```julia
readIJ(filename::AbstractString)
readIJ(filename::AbstractString, sep)
```

IO.jl:4



### readIJV
to read a simple edge list, each line being an (i, j, v) pair. The parens should not be there in the format, just commas separating. To generate this format in Matlab, you just need to be careful to write the vertex indices with sufficient precision.  For example, you can do this

```
>> [ai,aj,av] = find(triu(a));
>> dlmwrite('graph.txt',[ai,aj,av],'precision',9);
```


```julia
readIJV(filename::AbstractString)
```

IO.jl:25



### writeIJV
Writes the upper portion of a matrix in ijv format, one row for each edge, separated by commas.  Only writes the upper triangular portion. The result can be read from Matlab like this:

```
>> dl = dlmread('graph.txt');
>> a = sparse(dl(:,1),dl(:,2),dl(:,3));
>> n = max(size(a))
>> a(n,n) = 0;
>> a = a + a';
```


```julia
writeIJV(filename::AbstractString, mat)
```

IO.jl:52



