%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Marcus Tan on 6/8/2013
%%% Last modified date: 6/22/2014
%%% Copyright 2013 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This function numbers the edges of the triangular elements and gives the
% element number of the elements sharing each edge. 
function [edge_nodes, elem_edges] =find_edge_nodes_elem_edges(elem_nodes)
    edge_nodes = [elem_nodes(:,[1,2]); elem_nodes(:,[2,3]); elem_nodes(:,[3,1])];
    [edge_nodes,~,ic] = unique(sort(edge_nodes,2),'rows');
    nElems = size(elem_nodes,1);
    elem_edges = [ic(1:nElems),ic(nElems+1:2*nElems),ic(2*nElems+1:3*nElems)];   
end