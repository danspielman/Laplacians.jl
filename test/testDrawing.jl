if !haskey(ENV,"LAPLACIANS_NOPLOT")
    println("*** TESTING PLOTTING ***")
    a = wted_chimera(102,2)
    spectral_drawing(a)
else
    println("NOT TESTING PLOTTING")
end

