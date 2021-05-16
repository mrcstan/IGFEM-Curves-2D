%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Marcus Tan on 6/8/2013
%%% Last modified date: 6/8/2013
%%% Copyright 2013 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this function finds the edge number of the edges of each element
% Assume each row of edge_node has been sorted in ascending order
function elem_edges=find_elem_edge2(elem_nodes,edge_nodes)    
   elem_edges = int32(zeros(size(elem_nodes,1),3));
   [~,locb] = ismember(sort(elem_nodes(:,1:2),2),edge_nodes,'rows');  
   elem_edges(:,1) = locb;
   [~,locb] = ismember(sort(elem_nodes(:,2:3),2),edge_nodes,'rows');  
   elem_edges(:,2) = locb;
   [~,locb] = ismember(sort(elem_nodes(:,[1,3]),2),edge_nodes,'rows');  
   elem_edges(:,3) = locb;
end