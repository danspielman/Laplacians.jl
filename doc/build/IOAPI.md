
<a id='IO-1'></a>

# IO

<a id='Laplacians.readIJ' href='#Laplacians.readIJ'>#</a>
**`Laplacians.readIJ`** &mdash; *Function*.



To read a simple edge list, each line being an (i, j) pair


<a target='_blank' href='https://github.com/danspielman/Laplacians.jl/tree/ce307408e6b2290f943763f162deee084e2b9097/doc/../src/IO.jl#L2' class='documenter-source'>source</a><br>

<a id='Laplacians.readIJV-Tuple{AbstractString}' href='#Laplacians.readIJV-Tuple{AbstractString}'>#</a>
**`Laplacians.readIJV`** &mdash; *Method*.



To read a simple edge list, each line being an (i, j, v) pair. The parens should not be there in the format, just commas separating. To generate this format in Matlab, you just need to be careful to write the vertex indices with sufficient precision.  For example, you can do this

```
>> [ai,aj,av] = find(triu(a));
>> dlmwrite('graph.txt',[ai,aj,av],'precision',9);
```


<a target='_blank' href='https://github.com/danspielman/Laplacians.jl/tree/ce307408e6b2290f943763f162deee084e2b9097/doc/../src/IO.jl#L14-L24' class='documenter-source'>source</a><br>

<a id='Laplacians.writeIJV-Tuple{AbstractString,Any}' href='#Laplacians.writeIJV-Tuple{AbstractString,Any}'>#</a>
**`Laplacians.writeIJV`** &mdash; *Method*.



Writes the upper portion of a matrix in ijv format, one row for each edge, separated by commas.  Only writes the upper triangular portion. The result can be read from Matlab like this:

```
>> dl = dlmread('graph.txt');
>> a = sparse(dl(:,1),dl(:,2),dl(:,3));
>> n = max(size(a))
>> a(n,n) = 0;
>> a = a + a';
```


<a target='_blank' href='https://github.com/danspielman/Laplacians.jl/tree/ce307408e6b2290f943763f162deee084e2b9097/doc/../src/IO.jl#L38-L50' class='documenter-source'>source</a><br>

