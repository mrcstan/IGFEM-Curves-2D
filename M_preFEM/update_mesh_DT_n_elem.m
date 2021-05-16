%%% Created by Marcus Tan on 3/23/2016
%%% Copyright 2016 University of Illinois at Urbana-Champaign. 
%%% All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % This function constructs the Delaunay Triangulation object and
 % update other members of mesh.elem to be consistent with the Delaunay
 % Triangulation
function [DT,elem] = update_mesh_DT_n_elem(nodeCoords, ...
                                           edgeNodes, ...
                                           elem)
% reconstruct DT if refinement has been carried out
% or construct DT if DT is empty
DT = delaunayTriangulation(nodeCoords, ...
                           edgeNodes);

                       
if nargout > 1
    % update heat source and materials  
    [~,new2oldElnums] = ismember(sort(DT.ConnectivityList,2),...
                                 sort(elem.elem_node,2),'rows');
    elem.elem_node = DT.ConnectivityList;
    elem.dualedge = elems_sharing_edge_nodes(elem.elem_node,...
                                             size(nodeCoords,1)); 
    elem.heatSource = elem.heatSource(new2oldElnums);
    elem.material = elem.material(new2oldElnums);
end
end