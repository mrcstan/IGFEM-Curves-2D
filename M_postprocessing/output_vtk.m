outfile = '10x5.vtk';
UUR2 = update_enrichment_node_value(UUR,mesh.node,mesh.edge,mesh.elem);   
matlab2vtk_scalar(outfile,scalarname,mesh.elem,mesh.node,UUR2);