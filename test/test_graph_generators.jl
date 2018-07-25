

@test L.path_graph_ijv(1) == L.empty_graph_ijv(1)
@test path_graph(1) == empty_graph(1)
path_graph(5)


@test L.complete_graph_ijv(1) == L.empty_graph_ijv(1)
@test complete_graph(1) == empty_graph(1)
L.complete_graph_ijv(5)

@test L.ring_graph_ijv(1) == L.empty_graph_ijv(1)
@test ring_graph(1) == empty_graph(1)
ring_graph(5)

generalized_ring(11,[1 2 5])
@test generalized_ring(5,[5 10]) == empty_graph(5)
@test generalized_ring(1,[5 10]) == empty_graph(1)

rand_gen_ring(100,3)
rand_gen_ring(100,3,verbose=true)

@test hypercube(1) == path_graph(2)
hypercube(3)

a = L.complete_graph_ijv(5)
b = L.cbt_ijv(4)*3
ab = product_graph(a,b);

@test sparse(ab) == product_graph(sparse(a),sparse(b))

grid2(2,3,isotropy=2.0)
grid2(3)

grid3(3)
grid3(2,3,4)

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
L.ErdosRenyi_ijv(10,20)

srand(1)
L.ErdosRenyiCluster_ijv(6,30)
L.ErdosRenyiCluster_ijv(30,10)

srand(1)
L.ErdosRenyiClusterFix_ijv(20,20)
L.ErdosRenyiClusterFix_ijv(10,50)

L.firstn(L.grid2_ijv(10),45)


;




