scale.uo = 20;
scale.ud = 2000*0.1/(2*0.6);
maxTemp = (max(UUR2)-scale.uo)/scale.ud

aveTemp =  average_temp(mesh.elem,mesh.node.coords,UUR2);
aveTemp = (aveTemp-scale.uo)/scale.ud

aveEdgeLength = mean(mesh.edge.length)
minEdgeLength = min(mesh.edge.length)
maxEdgeLength = max(mesh.edge.length)

dof = mesh.node.n_node-mesh.node.Dirichlet.n_pre_temp