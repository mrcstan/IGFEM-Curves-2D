%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Marcus Tan on 6/10/2013
%%% Last modified date: 7/2/2014
%%% Copyright 2013 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This function set the Dirichlet and Neumann boundary conditions on a
% rectangular domain
% INPUT:
%   BCs.boundaries: an array of numel N of boundary numbers. 
%                   1: bottom, 2: right, 3: top, 4: left
%   BCs.types: a numel N array boundary condition type
%               1: Dirichlet, 2: Neumann
%   BCs.values_or_funcs: a numel Ncell array of values or function handles 
%                        or a mixed of both
% OUTPUT:
%   node.Dirichlet.n_pre_temp: total number of all the Dirichlet nodes
%   node.Dirichlet.temp_node:  all the Dirichlet nodes
%   node.Dirichlet.temp_value: the value at the Dirichlet nodes. 
%                            the value at the presribed enrichment nodes
%                            are set to zero. 
%   elem.Neumann.n_heatFlux
%   elem.Neumann.heatFlux_elem
%   elem.Neumann.heatFlux_value
%   elem.Neumann.heatFlux_surface
%  REMARK:
%  The offset value at an isolated Dirichlet enrichment node such as the inlet of
%  a channel cannot be defined as the values at the adjacent original nodes
%  are unknown
function [elem, node] = set_boundary_conditions(BCs,...
                                                elem,...
                                                node,...                                               
                                                boundary,...
                                                tol)
                                                                                        

nNodes = size(node.coords,1);
node.Dirichlet.n_pre_temp = 0;
node.Dirichlet.temp_node = zeros(nNodes,1);
node.Dirichlet.temp_value = zeros(nNodes,1);

nElems = size(elem.elem_node,1);
elem.Neumann.n_heatFlux = 0;
elem.Neumann.heatFlux_elem = zeros(nElems,1);
elem.Neumann.heatFlux_value = zeros(nElems,1);
elem.Neumann.heatFlux_surface = zeros(nElems,1);

auxSet1 = [1,1,2,2]; 
auxSet2 = [boundary.xi,boundary.xf,boundary.yi,boundary.yf];
auxSet3 = [1,2,3]; %originally S1,S2 and S3
nodes = (1:node.n_node)';

%{
nEdges = size(edge.edge_node,1);
d2p = sparse(edge.edge_node(:,[1,2]),...
             edge.edge_node(:,[2,1]),...
             [1:nEdges,1:nEdges],...
             node.nOriginalNode,node.nOriginalNode); 
%}  

for i = 1:numel(BCs.boundaries)
    bn = BCs.boundaries(i);
    bnNodes = nodes(abs(node.coords(1:node.nOriginalNode,auxSet1(bn)) - auxSet2(bn)) < tol);
     % Important to reset bnNodeVals !
    bnNodeVals = [];
    switch (BCs.types(i))
        case 1 % Dirichlet BC
            nPreTemp = node.Dirichlet.n_pre_temp + numel(bnNodes);
            node.Dirichlet.temp_node((node.Dirichlet.n_pre_temp+1):nPreTemp) ...
                    =  bnNodes;
            if (isa(BCs.values_or_funcs{i},'function_handle'))
                bnNodeVals(1:numel(bnNodes)) = BCs.values_or_funcs{i}(node.coords(bnNodes,1), ...
                                                    node.coords(bnNodes,2));
            else
                bnNodeVals(1:numel(bnNodes)) = BCs.values_or_funcs{i};
            end
            node.Dirichlet.temp_value((node.Dirichlet.n_pre_temp+1):nPreTemp) ...
                    =  bnNodeVals(:);
            node.Dirichlet.n_pre_temp = nPreTemp;
        case 2 % Neumann BC
            %{
            for j = 1:numel(bnNodes)
                if (bnNodes(i))
                    neighborNodes = find(d2p(bnNodes(i),:));
                    [neighborBnNode,locBnNodes] = ismember(neighborNodes,bnNodes);
                    bnNodes(locBnNodes) = 0;
                    neighborBnNode = neighborNodes(neighborBnNode)
                    edge = d2p(bnNodes(j),neighborBnNode);
                    
                end
            end
            %}
            elemsSharingNodePairs = elem.dualedge(bnNodes,bnNodes); % all elements that have the node pairs
            [row,col] = find(elemsSharingNodePairs); % find the node pairs that are connected by an edge
            for j = 1:numel(row)
                if (~elem.dualedge(bnNodes(col(j)),bnNodes(row(j))))
                    el = elem.dualedge(bnNodes(row(j)),bnNodes(col(j)));
                    elem.Neumann.n_heatFlux = elem.Neumann.n_heatFlux + 1;
                    elem.Neumann.heatFlux_elem(elem.Neumann.n_heatFlux) = el;
                    if (isa(BCs.values_or_funcs{i},'function_handle'))
                        Xc = mean(node.coords(bnNodes([row(j),col(j)]),:),1); % find midpoint
                        elem.Neumann.heatFlux_value(elem.Neumann.n_heatFlux) ...
                            = BCs.values_or_funcs{i}(Xc(1),Xc(2));
                    else
                        elem.Neumann.heatFlux_value(elem.Neumann.n_heatFlux) = BCs.values_or_funcs{i};
                    end
                    locEdge = local_edge(elem.elem_node(el,:),bnNodes([row(j),col(j)]));
                    elem.Neumann.heatFlux_surface(elem.Neumann.n_heatFlux) = auxSet3(locEdge);
                end
            end
    end
end
elem.Neumann.heatFlux_elem = elem.Neumann.heatFlux_elem(1:elem.Neumann.n_heatFlux);
elem.Neumann.heatFlux_value = elem.Neumann.heatFlux_value(1:elem.Neumann.n_heatFlux);
elem.Neumann.heatFlux_surface = elem.Neumann.heatFlux_surface(1:elem.Neumann.n_heatFlux);

% add constraint nodes as Dirichlet nodes if they are original nodes
del = false(numel(node.constraint.temp_node),1);
for i = 1:node.constraint.number
    if (node.constraint.temp_node(i) <= node.nOriginalNode)
        fprintf('\nTransferring constraint node %i to Dirichlet node\n', ...
                node.constraint.temp_node(i));    
        del(i) = true;
        node.Dirichlet.n_pre_temp = node.Dirichlet.n_pre_temp + 1;
        node.Dirichlet.temp_node(node.Dirichlet.n_pre_temp) = node.constraint.temp_node(i);
        node.Dirichlet.temp_value(node.Dirichlet.n_pre_temp) = node.constraint.temp_value(i);
    end
end
node.constraint.temp_node(del) = [];
node.constraint.temp_value(del) = [];
node.constraint.number = numel(node.constraint.temp_node);
node.Dirichlet.temp_node = node.Dirichlet.temp_node(1:node.Dirichlet.n_pre_temp);
node.Dirichlet.temp_value = node.Dirichlet.temp_value(1:node.Dirichlet.n_pre_temp);

% remove repeated Dirichlet nodes
if (node.Dirichlet.n_pre_temp > 2)
    [node.Dirichlet.temp_node,ia,~] = unique(node.Dirichlet.temp_node);
    if (numel(ia) ~= node.Dirichlet.n_pre_temp)

        warning(['existence of corner node where > 1 Dirichlet boundaries meets'...
                 ' or inlet temperature specified on Dirichlet boundary']) 
    end
    node.Dirichlet.temp_value = node.Dirichlet.temp_value(ia);
end
node.Dirichlet.n_pre_temp = numel(node.Dirichlet.temp_node);
end