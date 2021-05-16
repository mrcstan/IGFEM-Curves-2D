%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Marcus Tan on 10/4/2013
%%% Modified on 5/1/2014
%%% Copyright 2013 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this function flip an edge shared by two triangles
% WARNING:
%   Assume that check has already been done to make sure the quadrilateral 
%   formed by the two triangles is convex!
% INPUT: 
%   edge_node: length nx2 array of the nodes of the edges of the elements
%               where the ith row contains the nodes of the ith edge
%   elem_node: size 2x3 array of the connectivity of the 2 elements sharing
%              the edge
%   elem_edge (optional): size 2x3 array of the connectivity of the 2 elements sharing
%                         the edge
% OUTPUT:
%   edge_node: length 2 array of the nodes of the flipped edge
%   elem_node,elem_edge: updated connectivity
function [flip_edge_node,elem_node,elem_edge] = flip_edge(edge_node,elem_node,elem_edge)
      if (nargin > 2)
        [flip_edge,locFlipEdge1,locFlipEdge2] = intersect(elem_edge(1,:),elem_edge(2,:));
      else
          flip_edge = 1;
      end
      if(isempty(flip_edge))
          error('common edge not found')
      end
      node = unique(elem_node);      
      oppNode1 = setdiff(node,elem_node(1,:)); % node opposite to element 1
      oppNode2 = setdiff(node,elem_node(2,:)); % node opposite to element 2
      sharedLocNode1 = find(elem_node(1,:)==edge_node(flip_edge,1),1,'first'); % it is also possible to do find(elem_node(1,:)==edge_node(2))
      sharedLocNode2 = find(elem_node(2,:)==edge_node(flip_edge,2),1,'first');
      sharedNode1 = elem_node(1,sharedLocNode1);
      sharedNode2 = elem_node(2,sharedLocNode2);
      % update the nodes of the flipped edge
      flip_edge_node(1) = oppNode2;
      flip_edge_node(2) = oppNode1;
      % update connectivity of the nodes of the elements
      elem_node(1,sharedLocNode1) = oppNode1; % swapped the node shared by the edge with opposite node of elem 1
      elem_node(2,sharedLocNode2) = oppNode2; % swapped the node shared by the edge with opposite node of elem 2    
      
      if (nargin > 2)
          % update the connectivity of the edges of the elements      
          % find the edge opposite to sharedLocNode1
          oppLocEdge1 = ~(edge_node(elem_edge(1,:),1)== sharedNode2 | ...
                            edge_node(elem_edge(1,:),2)== sharedNode2); 
          oppEdge1 = elem_edge(1,oppLocEdge1);
          % find the edge opposite to sharedLocNode2
          oppLocEdge2 = ~(edge_node(elem_edge(2,:),1)== sharedNode1 |...
                            edge_node(elem_edge(2,:),2)== sharedNode1); 
          oppEdge2 = elem_edge(2,oppLocEdge2);
          elem_edge(1,locFlipEdge1) = oppEdge2;
          elem_edge(1,oppLocEdge1) = flip_edge;
          elem_edge(2,locFlipEdge2) = oppEdge1;
          elem_edge(2,oppLocEdge2) = flip_edge;
      end
end