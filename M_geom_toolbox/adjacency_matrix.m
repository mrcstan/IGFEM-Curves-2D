%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Marcus Tan on 3/3/2014
%%% Last modified date: 8/21/2014
%%% Copyright 2014 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this function computes the element to node, node to node, element to element adjacency matrices
% INPUT: 
%	elem_nodes: a nElems x nNodesPerElem array of the element connectivity, 
%               where nNodesPerElem is the number of nodes per element
%   nNodes: total number of nodes
% OUTPUT:
%	elem2node: nElems x nNodes array. elem2node(i,j) = 1 if elem i contains node j
%													 = 0 otherwise.
%	node2node: nNodes x nNodes array. node2node(i,j) = number of elements attached to node i when i=j,
%													 = number of elements attached to edge ij when nodes i and j are neighbours
%													 = 0 otherwise
%   elem2elem: nElems x nElems array. elem2elem(i,j) = number of nodes shared between elements i and j
function [elem2node,node2node,elem2elem] = adjacency_matrix(elem_nodes,nNodes)
nNodesPerElem = size(elem_nodes,2);
nElems = size(elem_nodes,1);
repeatElems = repmat(1:nElems,nNodesPerElem,1);
elem2node = sparse(repeatElems(:),elem_nodes',1,nElems,nNodes,nElems*nNodesPerElem);

if (nargout > 1)
    node2node = elem2node'*elem2node;
end
if (nargout > 2)
    elem2elem = elem2node*elem2node';
end
end