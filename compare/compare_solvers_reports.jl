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
  p = sortperm(dic["$(sol1)_$(str)"])
  for name in dic["names"]
      Plots.plot!(dic["$(name)_$(str)"][p], label=name)
  end
  display(plt)
end
