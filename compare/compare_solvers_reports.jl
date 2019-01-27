using Plots
pyplot()
#plotly()
#plotlyjs() #plotlyjs issue: https://github.com/JuliaPlots/Plots.jl/issues/1793
using Statistics
using Printf

function total_time(dic; names = dic["names"], labels = dic["names"])
    plt = Plots.plot(title = "Seconds Distribution")
    for i in 1:length(names)
        name = names[i]
        Plots.plot!(sort(dic["$(name)_tot"]), label=labels[i], linewidth=3)
    end
    display(plt)
    return plt
end

function total_time_relative(dic; names = dic["names"], labels=names, cols=distinguishable_colors(length(names)))

    if isa(cols,Dict)
        cols = [cols[x] for x in names];
    end

  plt = Plots.plot(title = "(Time / Time for $(labels[1])), Distribution")
  rel = dic["$(names[1])_tot"]

  mn = 1
  mx = 1
    
  for i in 2:length(names)
        name = names[i]
      y = sort(min.(dic["$(name)_tot"] ./ rel, 10))
      mn = min(mn,minimum(y))
      mx = max(mx,maximum(y))
      Plots.plot!(y, label=labels[i], linewidth=3, color=cols[i])
  end
  Plots.plot!(ylimits=(0:10),yticks=[0;round(mn,2);1:10])

  plot!(ytickfont = font(16))
  plot!(legendfont = font(16))
  plot!(guidefont = font(16))
    plot!(xticks=[])


  display(plt)
  return plt
end

function total_time_relative_unsort(dic; names = dic["names"], labels=names, cols=distinguishable_colors(length(names)))

    if isa(cols,Dict)
        cols = [cols[x] for x in names];
    end

  plt = Plots.plot(title = "(Time / Time for $(labels[1])), Distribution")
  rel = dic["$(names[1])_tot"]

  mn = 1
  mx = 1
    
  for i in 2:length(names)
        name = names[i]
      #y = sort(min.(dic["$(name)_tot"] ./ rel, 10))
      y = min.(dic["$(name)_tot"] ./ rel, 10)
      mn = min(mn,minimum(y))
      mx = max(mx,maximum(y))
      Plots.plot!(y, label=labels[i], linewidth=3, color=cols[i])
  end
  Plots.plot!(ylimits=(0:10),yticks=[0;round(mn,2);1:10])

  plot!(ytickfont = font(16))
  plot!(legendfont = font(16))
  plot!(guidefont = font(16))
    plot!(xticks=[])


  display(plt)
  return plt
end



function total_time_nnz(dic; names = dic["names"],
                        labels = names,
                        cols=distinguishable_colors(length(names)),
                        ref=[])
    #dic = copy(dic)
    if isa(cols,Dict)
       cols = [cols[x] for x in names];
    end

    if ~isempty(ref)
        yref = 10*dic["$(ref)_tot"] ./ dic["ne"]
    end
    
    plt = Plots.plot()
    mn = Inf
    mx = 0
    for name in names
        y = dic["$(name)_tot"] ./ dic["ne"]
        if ~isempty(ref)
            y = min.(y, yref)
        end
            
        mn = min(mn,minimum(y))
        mx = max(mx,maximum(y[y .< Inf]))
    end
    b = floor(log(mx)/log(10))
    for i in 1:length(names)
        name = names[i]
               y = dic["$(name)_tot"] ./ dic["ne"]
        if ~isempty(ref)
            y = min.(y, yref)
        end
      y = sort(y)
      y[y .== Inf] .= mx
      y = y * 10^(-b)  
      Plots.plot!(y, label=labels[i], linewidth=3, color=cols[i])
    end
  plot!(ytickfont = font(12))
  plot!(legendfont = font(16))
  plot!(guidefont = font(16))
    plot!(xticks=[])
  title!("$(10.0^(round.(Int,b))) Seconds / Edge, Distribution")
  display(plt)
  return plt, mn, mx
end

function total_time_nnz_loglog(dic; names = dic["names"],
                        labels = names,
                        cols=distinguishable_colors(length(names)),
                            #ref=[]
                            )
    if isa(cols,Dict)
       cols = [cols[x] for x in names];
    end

    # if ~isempty(ref)
    #     yref = 10*dic["$(ref)_tot"] ./ dic["ne"]
    # end

    x = dic["ne"]
    
    plt = Plots.plot()
    mn = Inf
    mx = 0
    for name in names
        y = dic["$(name)_tot"]
        # if ~isempty(ref)
        #     y = min.(y, yref)
        # end
        mn = min(mn,minimum(y))
        if any(y .< Inf)
            mx = max(mx,maximum(y[y .< Inf]))
        end
    end

   # b = floor(log(mx)/log(10))
    for i in 1:length(names)
        name = names[i]
        y = dic["$(name)_tot"]
        # if ~isempty(ref)
        #     y = min.(y, yref)
        # end
        #y = sort(y)
        y[y .== Inf] .= 2*mx
        #y = y * 10^(-b)  
        Plots.plot!(x, y, label="", linewidth=3, color=cols[i])
        Plots.scatter!(x, y, label=labels[i], linewidth=3, color=cols[i])
    end
  Plots.plot!(xaxis=(:log), yaxis=(:log))
  plot!(ytickfont = font(12))
  plot!(legendfont = font(16))
    plot!(guidefont = font(16))
    #plot!(xticks=[])
  #title!("$(10.0^(round.(Int,b))) Seconds / Edge, Distribution")
  display(plt)
  return plt, mn, mx
end

function total_pair(dic, name1, name2)
    plt = Plots.plot()
    mn = Inf
    mx = 0
    y1 = dic["$(name1)_tot"]
    y2 = dic["$(name2)_tot"]
    
    mn = min(mn,minimum(y1),minimum(y2))
    if any(y1 .< Inf)
        mx = max(mx,maximum(y1))
    end
    if any(y2 .< Inf)
        mx = max(mx,maximum(y2))
    end
    y1[y1 .== Inf] .= 2*mx
    y2[y2 .== Inf] .= 2*mx

    Plots.plot!([mn,mx], [mn,mx], labels="", linewidth=3)
    Plots.scatter!(y1, y2, label="", linewidth=3)

    
    Plots.plot!(xaxis=(:log), yaxis=(:log))
    plot!(ytickfont = font(12))
    plot!(legendfont = font(16))
    plot!(guidefont = font(16))
    display(plt)
    return plt, mn, mx
end






function time_vs_n(dic; names = dic["names"], labels=names, cols=distinguishable_colors(length(names)))
    if isa(cols,Dict)
        cols = [cols[x] for x in names];
    end
    
    plt = Plots.plot()
    for i in 1:length(names)
        name = names[i]
      Plots.plot!(dic["nv"],dic["$(name)_tot"], label=labels[i], linewidth=3, color=cols[i])
  end
    Plots.plot!(xaxis=(:log), yaxis=(:log))
    Plots.plot!(xlabel="Number Vertices",ylabel="Seconds")
      plot!(ytickfont = font(12))
      plot!(xtickfont = font(12))
  plot!(legendfont = font(16))
  plot!(guidefont = font(16))
  display(plt)
  return plt
end


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

  p1 = plot_relative(dic,"tot","Total Time")
  p2 = plot_relative(dic,"build","Build Time")
  p3 = plot_relative(dic,"solve","Solve Time")
  p4 = plot_relative(dic,"its","Iterations")
  return plot(p1,p2,p3,p4)
end

function plot_relative_unsort(dic)

  p1 = plot_relative_unsort(dic,"tot","Total Time")
  p2 = plot_relative_unsort(dic,"build","Build Time")
  p3 = plot_relative_unsort(dic,"solve","Solve Time")
  p4 = plot_relative_unsort(dic,"its","Iterations")
  return plot(p1,p2,p3,p4)
end

function vec_stats(v)
    v = sort(vec(v))
    n = length(v)
    @printf("mean: %e \n",mean(v))
    pts = round.(Int, [1, n/4, n/2, 3n/4, n],RoundUp)
    @printf(" quarts: %e %e %e %e %e ", v[pts[1]], v[pts[2]], v[pts[3]], v[pts[4]], v[pts[5]])
    println()
end

function vec_stats(v,cnt)
    v = sort(vec(v))
    n = length(v)
    @printf("mean: %e \n",mean(v))
    pts = round.(Int, [1, n/4, n/2, 3n/4, n],RoundUp)
    @printf(" quarts: %e %e %e %e %e ", v[pts[1]], v[pts[2]], v[pts[3]], v[pts[4]], v[pts[5]])
    print("  (Inf: ",sum(isinf.(cnt)), ")")

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
        threshrow = min.(thisrow, 10*baserow)
        Plots.plot!(threshrow, label=name)

        print(key, " ")
        vec_stats(threshrow, thisrow)
  end
  display(plt)
end



function plot_relative(dic,str,title)

    plt = scatter(title = title)
    sol1 = dic["names"][1]
    key = "$(sol1)_$(str)"

    baserow = dic["$(sol1)_$(str)"]

    for i in 2:length(dic["names"])
        name = dic["names"][i]
        rowname = "$name / $(sol1)"
        key = "$(name)_$(str)"
        thisrow = dic[key]
        threshrow = min.(thisrow, 10*baserow) ./ baserow
        scatter!(sort(threshrow), label=rowname)

        print(str, " " ,rowname, " ")
        vec_stats(threshrow, thisrow)
    end
    #display(plt)
    return plt
end

function plot_relative_unsort(dic,str,title)

    plt = scatter(title = title)
    sol1 = dic["names"][1]
    key = "$(sol1)_$(str)"

    baserow = dic["$(sol1)_$(str)"]

    for i in 2:length(dic["names"])
        name = dic["names"][i]
        rowname = "$name / $(sol1)"
        key = "$(name)_$(str)"
        thisrow = dic[key]
        threshrow = min.(thisrow, 10*baserow) ./ baserow
        scatter!(sort(threshrow), label=rowname)

        print(str, " " ,rowname, " ")
        vec_stats(threshrow, thisrow)
    end
    #display(plt)
    return plt
end

function compare_two(dic,sol1,sol2)

  plt = Plots.plot(title = "Total Time, Relative")
    key1 = "$(sol1)_tot"
    key2 = "$(sol2)_tot"

    tot1 = dic[key1]
    tot2 = dic[key2]

    Plots.plot!(sort(tot1./tot2), label="$(sol1) / $(sol2)")
    Plots.plot!(sort(tot2./tot1), label="$(sol2) / $(sol1)")

  display(plt)

  plt = Plots.plot(title = "Total Time / nnz")
    key1 = "$(sol1)_tot"
    key2 = "$(sol2)_tot"

    tot1 = dic[key1]
    tot2 = dic[key2]
    ne = dic["ne"]

    Plots.plot!(sort(tot1./ne), label="$(sol1) / nnz")
    Plots.plot!(sort(tot2./ne), label="$(sol2) / nnz")

  display(plt)
end

"""
    sd = mask_dict(mcf,mask)

For example:
```
mask = mcf["nv"] .> 20000
@show sum(mask)
sd = mask_dict(mcf,mask)
```
"""
function mask_dict(dic,mask)
    sol = dic["names"][1]
    ntests = length(dic["$(sol)_tot"])    
    
    subdict = Dict()

    for (key, value) in dic
        if length(dic[key]) == ntests
            subdict[key] = dic[key][mask]
        else
            subdict[key] = dic[key]
        end
    end
    return subdict
end
