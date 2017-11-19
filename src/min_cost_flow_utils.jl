
#=   
Find the the inactive edges based on the complementarity conditions. The exact procedure is 
descriped in the pdnet paper.
=#
# function find_inactive_edges(x,s,z,zeta0)
# m = size(x)[1]
# xT = zeros(size(x))
# sT = zeros(size(s))
# zT = zeros(size(z))
# inact_edges = zeros(m)

# for i in 1:m

#   if x[i]/s[i] < zeta0  && (u[i]-x[i])/z[i] > 1/zeta0
#     xT[i] = 0
#     inact_edges[i] = -1
#   elseif x[i]/s[i] > 1/zeta0 && (u[i]-x[i])/z[i] < zeta0
#     xT[i] = u[i]
#     inact_edges[i] = 1

# end

# zeta0 = zeta0 * 0.95

# return (xT,sT,zT,zeta0,inact_edges)

# end

#=   
Computes the maximum weight forest restricted to inact_edges. This implements Kruskal's algoirthm.
=#

function max_wt_tree_opt_face(mat,d,inact_edges)
    n = size(mat)[1] 
  (ai,aj,av) = findnz(mat)
  comps = IntDisjointSets(n)
  ord = sortperm(d, rev=true)
  treeinds = zeros(Int,n-1)
  numintree = 0

  for i in ord
        if  (abs(inact_edges[i])< 0.5) && !DataStructures.in_same_set(comps,ai[i],aj[i])
      numintree = numintree+1
      treeinds[numintree] = i
      DataStructures.union!(comps,ai[i],aj[i])
    end
  end

  vert_ind = zeros(n)
  for i in 1:n
    vert_ind[i]=find_root(comps,i)
  end

  comp_ind = unique(vert_ind)
  treeinds = treeinds[1:end-length(comp_ind)+1]
  forest = sparse(ai[treeinds],aj[treeinds],av[treeinds],n,n)
  forest = tree + tree'

    return (comp_ind, vert_ind,forest)
end


#=   
This projects the dual variable y onto a dual face defined by a forest.
Given a forest F, this computes arg_min_y ||y - y0||^2  s.t. B_F*y = c_F.
=#
# function proj_y_opt_face(inactive)

# end

# function y_opt_face()

# end

# function build_max_flow_instance(B,b,xu,xl)

#   bnew = b - B*xu
#   #dadf

# end

#= Checks for feasibiligy of B^Tx = b, x>=0, u>=x.
We do this using maxflow
=#

function check_flow_feasibility(edges,caps,dems)

n = size(dems,1)
m = size(edges,1)
anew = spzeros(n+2,n+2)
demnew = zeros(n+2)


demnew[3:end] = dems

dem_pos = find(dems .> 0)
n_pos = size(dem_pos,1) 
dem_neg = find(dems .< 0)
n_neg = size(dem_neg,1)

edgesnew = zeros(Int,m+n_pos+n_neg,2)
capsnew = zeros(m+n_pos+n_neg)
edgesnew[1:m,:] = edges + 2
capsnew[1:m] = caps

edgesnew[m+1:m+n_pos,:] = hcat(ones(n_pos), 2+dem_pos)
capsnew[m+1:m+n_pos]  = dems[dem_pos]

edgesnew[m+n_pos+1:end,:] = hcat(2+dem_neg, 2*ones(n_neg))
capsnew[m+n_pos+1:end] = -dems[dem_neg]



anew = sparse(edgesnew[:,1],edgesnew[:,2],capsnew,n+2,n+2)
flowGraph = DiGraph(anew)
F, f = maximum_flow(flowGraph, 1, 2,anew,algorithm=EdmondsKarpAlgorithm())

if F < sum(dems[dem_pos])
  return (false,F,f)
else 
  return (true,F,f) 

end



end
