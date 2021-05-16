%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Marcus Tan on 6/8/2013
%%% Last modified date: 6/8/2013
%%% Copyright 2013 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this function finds the edge number of the edges of each element
% Assume each row of edge_node has been sorted in ascending order
function elem_edges=find_elem_edge(elem_nodes,edge_nodes,nNodes)  
   nEdges = size(edge_nodes,1);
   d2p = sparse(edge_nodes(:,[1,2]),...
                edge_nodes(:,[2,1]),...
                [1:nEdges,1:nEdges],...
                nNodes,nNodes); 
   elem_edges = zeros(size(elem_nodes,1),3);
   ind = sub2ind(size(d2p),elem_nodes(:,1),elem_nodes(:,2));
   elem_edges(:,1) = d2p(ind);
   ind = sub2ind(size(d2p),elem_nodes(:,2),elem_nodes(:,3));
   elem_edges(:,2) = d2p(ind);
   ind = sub2ind(size(d2p),elem_nodes(:,1),elem_nodes(:,3));
   elem_edges(:,3) = d2p(ind);
end