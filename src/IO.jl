
"""to read a simple edge list, each line being an (i, j) pair"""
function readIJ(filename::AbstractString)
  edlist = readdlm(filename,',')
  n = maximum(edlist)
  m = size(edlist)
  edlist = convert(Array{Int64,2}, edlist)
  a = sparse(edlist[:,1],edlist[:,2],ones(m[1]),n,n)
  a = a + a'

  return a
end # readIJ

"""to read a simple edge list, each line being an (i, j, v) pair.
The parens should not be there in the format, just commas separating.
To generate this format in Matlab, you just need to be careful to write the
vertex indices with sufficient precision.  For example, you can do this

```
>> [ai,aj,av] = find(triu(a));
>> dlmwrite('graph.txt',[ai,aj,av],'precision',9);
```
"""
function readIJV(filename::AbstractString)
  data = readdlm(filename,',')
  n = maximum(data[:,1:2])
  m = size(data)
  edlist = convert(Array{Int64,2}, data[:,1:2])
  wts = convert(Array{Float64,1}, data[:,3])

  a = sparse(edlist[:,1],edlist[:,2],wts,n,n)
  a = a + a'

  return a
end # readIJV


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
