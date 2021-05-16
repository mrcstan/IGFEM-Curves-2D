%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Marcus Tan on 6/8/2013
%%% Last modified date: 6/8/2013
%%% Copyright 2013 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This function numbers the edges of the triangular elements and gives the
% element number of the elements sharing each edge. 
function edge_node=find_edge_node1(elem_node,nodeCoord)
    TR = TriRep(double(elem_node), nodeCoord(:,1), nodeCoord(:,2));
    edge_node=int32(edges(TR));
end