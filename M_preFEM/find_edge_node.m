%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Marcus Tan on 6/8/2013
%%% Last modified date: 6/22/2014
%%% Copyright 2013 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This function numbers the edges of the triangular elements and gives the
% element number of the elements sharing each edge. 
function edge_node=find_edge_node(elem_node)
    edge_node = [elem_node(:,[1,2]); elem_node(:,[2,3]); elem_node(:,[3,1])];
    edge_node = unique(sort(edge_node,2),'rows');
end