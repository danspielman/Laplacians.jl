#using Laplacians


function buildAPI(fileName, mod)
    noDocString = "No documentation found."

    v = names(mod)


    fh = open(fileName,"w")

    x = 0

    for i in 1:length(v)
        sym = v[i]
        x = eval(mod,sym)
        println(string(x))
        println(fh, "### ", string(sym))
#=        if (isa(x,Function))
          println(string(functionloc(x)))
        end
=#
       docmd = @doc(x)

        docstr = stringmime("text/plain", docmd )
        if (length(docstr) < 23) || (docstr[1:23] != noDocString)
            println(fh, docstr)
        else
            println("no docs for : ", sym)
        end

        extraInfo(fh, x)
        println(fh, "\n")

        
    end
    close(fh)
end


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







