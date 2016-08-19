using Laplacians


fileName = "docs/API/wholeAPI.md"

function extraInfo(x)
    fh = IOBuffer()
    if isa(x,Function)
        mt = methods(x)
        println(fh,"\n```julia")
        loc = " "
        firstit = true
        for meth in mt
            str = string(meth)
            ind = rsearchindex(str," at ")
            println(fh, str[1:(ind-1)])
            if firstit
                ind2 = rsearchindex(str,"/")+1
                loc = str[ind2:end]
                firstit = false
            end
        end
        println(fh,"```\n")

        println(fh, loc)

    elseif isa(x,DataType)
        str = Docs.typesummary(x)
        writemime(fh, "text/plain", str)
    end

    str = takebuf_string(fh)
end




noDocString = "No documentation found."

nms = names(Laplacians)

n = length(nms)

docstrs = Array(AbstractString,n)
linenums = zeros(Int64,n)
fileins = Array(AbstractString,n)
extras = Array(AbstractString,n)
strings = Array(AbstractString,n)

for i in 1:length(nms)
    sym = nms[i]
    x = eval(sym)

    strings[i] = string(sym)

    try
    
        if isa(x,Function)
            mt = methods(x)
            file, line = functionloc(mt.defs)
            fileins[i] = rsplit(file,"/")[end]
            linenums[i] = line
        end

        docmd = @doc string(x)
        
        docstr = stringmime("text/plain", docmd )
        if (length(docstr) < 23) || (docstr[1:23] != noDocString)
            docstrs[i] = docstr;
        end

        extras[i] = extraInfo(x)

    catch
        println("didn't generate docs for ", string(sym))
    end
        
end


ind = zeros(Bool,length(fileins))
for i in 1:length(fileins)
    ind[i] = isdefined(fileins,i)
end

u = unique(fileins[ind]);

fileins[!ind] = ""

for j in 1:length(u)
    infile = find(fileins .== u[j])
    infile = infile[sortperm(linenums[infile])]

    fn = split(u[j],'.')[1]
    
    fileName = "docs/API/" * fn * "API.md"
    fh = open(fileName, "w")

    println(fh, "# ", fn)
    
    for i in infile
        if isdefined(docstrs,i)
            println(fh, "### ", strings[i])
            println(fh, docstrs[i])
            println(fh, extras[i])
            println(fh,"\n")
        end
    end

    close(fh)

end


