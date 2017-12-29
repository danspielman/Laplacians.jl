if !haskey(ENV,"LAPLACIANS_NOPLOT")
    println("*** TESTING PLOTTING ***")
    a = wtedChimera(102,2)
    spectralDrawing(a)
else
    println("NOT TESTING PLOTTING")
end

