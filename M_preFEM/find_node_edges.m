%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Marcus Tan on 1/5/2014
%%% Last modified date: 3/15/2014
%%% Copyright 2014 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this function finds the edges sharing the nodes
function node_edges = find_node_edges(edge_nodes,nNodes)

node_edges = cell(nNodes,1);

for i = 1:nNodes
   [row,~] = find(edge_nodes == i);
  
   node_edges{i} = int32(row); 
end
end