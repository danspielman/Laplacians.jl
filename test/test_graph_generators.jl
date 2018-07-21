
@assert path_graph_ijv(1) == empty_graph_ijv(1)
@assert path_graph(1) == empty_graph(1)
path_graph(5)


@assert complete_graph_ijv(1) == empty_graph_ijv(1)
@assert complete_graph(1) == empty_graph(1)
complete_graph_ijv(5)

@assert ring_graph_ijv(1) == empty_graph_ijv(1)
@assert ring_graph(1) == empty_graph(1)
ring_graph(5)

generalized_ring(11,[1 2 5])
@assert generalized_ring(5,[5 10]) == empty_graph(5)
@assert generalized_ring(1,[5 10]) == empty_graph(1)

rand_gen_ring(100,3)
rand_gen_ring(100,3,verbose=true)

@assert hypercube(1) == path_graph(2)
hypercube(3)

a = complete_graph_ijv(5)
b = cbt_ijv(4)*3
ab = product_graph(a,b);

@assert sparse(ab) == productGraph(sparse(a),sparse(b))

grid2(2,3,isotropy=2.0)
grid2(3)

grid3(3)
grid3(2,3,4)

wgrid2(3)

srand(1)
rand_matching(4)
rand_matching(5)

srand(1)
rand_regular(7,1)
rand_regular(7,2)
rand_regular(6,1)
rand_regular(6,3)

srand(1)
grown_graph(11,1)
grown_graph(12,2)

srand(1)
grown_graph_d(50,2)
grown_graph_d(10,8)

srand(1)
pref_attach(6,5,0.5)
randperm(pref_attach(30,2,0.3))

srand(1)
ErdosRenyi_ijv(10,20)

srand(1)
ErdosRenyiCluster_ijv(6,30)
ErdosRenyiCluster_ijv(30,10)

srand(1)
ErdosRenyiClusterFix_ijv(20,20)
ErdosRenyiClusterFix_ijv(10,50)

firstn(grid2_ijv(10),45)


;




