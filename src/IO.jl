"""
    adj = read_graph(fn)

Read a graph from a file in IJ or IJV format.
That is, each line of the file should represent one edge.
Each edge should be specified by the indices of its vertices,
separated by a whitespace or a comma.  If the graph is weighted,
the weight should follow the second index.  For example, the unweighted
complete graph on 3 vertices would appear as the file

```
1 2
1 3
2 3
```

A weighted path on 3 vertices with edge weights `1.5` and `2.5` would be

```
1, 2, 1.5
2, 3, 2.5
```

The function tries to detect the delimiter type (comma or whitespace)
from the first line of the file.  The format must be consistent.
Vertex indices start at `1`.
"""
function read_graph(fn::AbstractString)

    # first, detect the delimiter type
    fh = open(fn)
    ln = readline(fh)
    close(fh)

    if occursin(",", ln)
        data = readdlm(fn,',')
    else
        data = readdlm(fn)
    end



    n = maximum(data[:,1:2])

    edlist = convert(Array{Int64,2}, data[:,1:2])

    if size(data,2) == 3
        wts = convert(Array{Float64,1},data[:,3])
    else
        wts = ones(size(data,1))
    end

    a = sparse(edlist[:,1],edlist[:,2],wts,n,n)
    a = a + a'

    return a
end




"""Writes the upper portion of a matrix in ijv format, one row for each edge,
separated by commas.  Only writes the upper triangular portion.
The result can be read from Matlab like this:

```
>> dl = dlmread('graph.txt');
>> a = sparse(dl(:,1),dl(:,2),dl(:,3));
>> n = max(size(a))
>> a(n,n) = 0;
>> a = a + a';
```
"""
function writeIJV(filename::AbstractString, mat)

  (ai,aj,av) = findnz(triu(mat))
  fh = open(filename,"w")
  for i in 1:length(ai)
    write(fh, "$(ai[i]),$(aj[i]),$(av[i])\n")
  end
  close(fh)

end #write IJV
