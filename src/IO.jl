
"""To read a simple edge list, each line being an (i, j) pair"""
function readIJ(filename::AbstractString, sep=',')
  edlist = readdlm(filename, sep)
  n = maximum(edlist)
  m = size(edlist)
  edlist = convert(Array{Int64,2}, edlist)
  a = sparse(edlist[:,1],edlist[:,2],ones(m[1]),n,n)
  a = a + a'

  return a
end # readIJ

"""To read a simple edge list, each line being an (i, j, v) pair.
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




" Writes a sparse matrix to a given text file "
function writeToFile{Tv,Ti}(filename::ASCIIString, a::SparseMatrixCSC{Tv,Ti})
    f = open(filename, "w")
    
    println(f, a.n, " ", a.m, " ", length(a.nzval))

    mat = tril(a);

    pos = 1
    for i in 1:length(mat.nzval)
        while mat.colptr[pos + 1] <= i
            pos = pos + 1
        end
        
        println(f, mat.rowval[i], " ", pos, " ", mat.nzval[i])
    end
    
    close(f)
end

" Reads a spare matrix from a given text file "
function readFromFile(filename::ASCIIString)
    r = readdlm(filename, ' ')
    
    m::Int64 = r[1,1]
    n::Int64 = r[1,2]
    nz::Int64 = r[1,3]
    U = Int64[]
    V = Int64[]
    W = Float64[]
    
    for i in 2:(nz+1)
        push!(U, r[i,1])
        push!(V, r[i,2])
        push!(W, r[i,3])

        push!(U, r[i,2])
        push!(V, r[i,1])
        push!(W, r[i,3])
    end
    
    return sparse(U,V,W,m,n)
end


