clear
close all
nx = 3;
ny = 3;
xi = 0.0;
xf = 1.0;
yi = 0.0;
yf = 1.0;
[elem_nodes,coords] = generate_uniform_mesh([nx,ny],[xi,xf,yi,yf],2);
%trimesh(elem_nodes,coords(:,1),coords(:,2),'color','k');
plot_mesh_labels(coords,elem_nodes,true,true)