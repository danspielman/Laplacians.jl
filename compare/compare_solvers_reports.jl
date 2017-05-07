using Plots
plotlyjs()

"""
    Given the list of solvers that were used, and a dic containing the results.
    Plot some graphs explaining how it went.
"""
function plots_and_stats(solvers,dic)

  plots_and_stats_raw(solvers,dic,"tot","Total Time")
  plots_and_stats_raw(solvers,dic,"build","Build Time")
  plots_and_stats_raw(solvers,dic,"solve","Solve Time")
  plots_and_stats_raw(solvers,dic,"its","Iterations")

end


function plots_and_stats(dic)

  plots_and_stats_raw(dic,"tot","Total Time")
  plots_and_stats_raw(dic,"build","Build Time")
  plots_and_stats_raw(dic,"solve","Solve Time")
  plots_and_stats_raw(dic,"its","Iterations")

end

function plot_relative(dic)

  plot_relative(dic,"tot","Total Time")
  plot_relative(dic,"build","Build Time")
  plot_relative(dic,"solve","Solve Time")
  plot_relative(dic,"its","Iterations")

end

function vec_stats(v)
    v = sort(v)
    n = length(v)
    print("mean: ",round(mean(v),3))
    pts = round(Int, [1, n/4, n/2, 3n/4, n])
    print(" quarts: ", v[pts])
    println()
end

function vec_stats(v,cnt)
    v = sort(v)
    n = length(v)
    print("mean: ",round(mean(v),3))
    pts = round(Int, [1, n/4, n/2, 3n/4, n])
    print(" quarts: ", v[pts])

    print("  (Inf: ",sum(isinf(cnt)), ")")

    println()
end

function plots_and_stats_raw(solvers,dic,str,title)

  plt = Plots.plot(title = title)
  sol1 = solvers[1]
  p = sortperm(dic["$(sol1.name)_$(str)"])
  for solver in solvers
      Plots.plot!(dic["$(solver.name)_$(str)"][p], label=solver.name)
  end
  display(plt)
end

function plots_and_stats_raw(dic,str,title)

  plt = Plots.plot(title = title)
    sol1 = dic["names"][1]
    key = "$(sol1)_$(str)"
    p = sortperm(dic[key])

    baserow = dic["$(sol1)_$(str)"][p]
    
    for name in dic["names"]
        key = "$(name)_$(str)"
        thisrow = dic[key][p]
        threshrow = min(thisrow, 10*baserow)
        Plots.plot!(threshrow, label=name)

        print(key, " ")
        vec_stats(threshrow, thisrow)
  end
  display(plt)
end



function plot_relative(dic,str,title)

  plt = Plots.plot(title = title)
    sol1 = dic["names"][1]
    key = "$(sol1)_$(str)"

    baserow = dic["$(sol1)_$(str)"]

    for i in 2:length(dic["names"])
        name = dic["names"][i]
        rowname = "$name / $(sol1)"
        key = "$(name)_$(str)"
        thisrow = dic[key]
        threshrow = min(thisrow, 10*baserow) ./ baserow
        Plots.plot!(sort(threshrow), label=rowname)

        print(str, " " ,rowname, " ")
        vec_stats(threshrow, thisrow)
  end
  display(plt)
end
