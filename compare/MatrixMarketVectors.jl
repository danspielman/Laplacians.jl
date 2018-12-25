struct ParseError
    error :: String
end

_parseint(x) = parse(Int, x)

function mmwritevec(filename, vec :: Array{Float64,1})
  open(filename, "w") do file
      # write mm header
      write(file, "%%MatrixMarket matrix array real general\n")

      # write matrix size and number of nonzeros
      write(file, "$(size(vec, 1)) $(size(vec, 2))\n")

      for i in 1:(size(vec, 1)-1)
          write(file, "$(vec[i])")
          write(file, "\n")
      end
      write(file, "$(vec[size(vec, 1)])")
  end
end

function mmreadvec(filename)
    open(filename,"r") do mmfile
        # Read first line
        firstline = chomp(readline(mmfile))
        tokens = split(firstline)
        if length(tokens) != 5
            throw(ParseError(string("Not enough words on first line: ", firstline)))
        end
        if tokens[1] != "%%MatrixMarket"
            throw(ParseError(string("Not a valid MatrixMarket header:", firstline)))
        end
        (head1, rep, field, symm) = map(lowercase, tokens[2:5])
        if head1 != "matrix"
            throw(ParseError("Unknown MatrixMarket data type: $head1 (only \"matrix\" is supported)"))
        end
        if !(rep == "array" && field == "real" && symm == "general")
            throw(ParseError("only   >array real general<  supported for now"))
        end
        
        # Skip all comments and empty lines
        ll   = readline(mmfile)
        while length(chomp(ll))==0 || (length(ll) > 0 && ll[1] == '%')
            ll = readline(mmfile)
        end
        # Read matrix dimensions (and number of entries) from first non-comment line
        dd = map(_parseint, split(ll))
        if length(dd) != 2 || dd[2] != 1
            throw(ParseError(string("Dim must be n by 1. Could not read in matrix dimensions from line:", ll)))
        end
        entries = dd[1]
        xx = Vector{Float64}(undef, entries)
        for i in 1:entries
            line = readline(mmfile)
            xx[i] = parse(Float64, line)
            # unsure if this breaks if last line empty
        end
        return xx
    end
end
