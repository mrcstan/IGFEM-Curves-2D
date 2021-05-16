path(path,'../M_preFEM')
tol = 1e-8;
for i = 1:size(mesh.node.coords,1)
    [ia,ib] = ismemberf(mesh.node.coords(i,:), ...
                        mesh.node.coords(i+1:end,:),'rows',tol);
    if ia
       error('node %i coincides with nodes %i',ia,ib)
    end
    
end