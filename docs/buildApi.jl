using Laplacians


fileName = "docs/API/wholeAPI.md"

function extraInfo(fh, x)
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
                loc = str[ind:end]
                firstit = false
            end
        end
        println(fh,"```\n")

        println(fh, loc)

    elseif isa(x,DataType)
        str = Docs.typesummary(x)
        writemime(fh, "text/plain", str)
    end
end




noDocString = "No documentation found."

v = names(Laplacians)


fh = open(fileName,"w")

x = eval(v[1])

for i in 1:length(v)
    sym = v[i]
    x = eval(sym)

    println(fh, "### ", string(sym))

    try
    
        docmd = @doc(x)
        
        docstr = stringmime("text/plain", docmd )
        if (length(docstr) < 23) || (docstr[1:23] != noDocString)
            println(fh, docstr)
        else
            println("no docs for : ", sym)
        end

        extraInfo(fh, x)

    catch
    end
        
    println(fh, "\n")


        
end
close(fh)




