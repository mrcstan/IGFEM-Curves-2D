% Taken from the paper "Short Bisection Implementation in MATLAB" by Long
% Chen
% Modified by Marcus Tan on June 22, 2014
%   INPUT:
%       elemMarkers: a length nElem logical array that specifies whether
%       dualedge: a nNodes x nNodes sparse matrix where each (i,j) entry is
%                 the element with an edge from node i to node j CCW. 
%                 The (j,i) entry is the neighbor element sharing the edge.
%                 This matrix is not symmetric.
%       d2p: a nNode x nNodes sparse matrix where each (i,j) entry is the
%            edge number with nodes i and j. This matrix is symmetric.
%   OUTPUT:
%       modifiedEdge: a length nEdge logical that indicates whether an edge
%                     is new or has been modified
%       elemProperty: a cell array, each cell is a length nElem array with
%                     element properties
function [node_coords,elem_nodes,edge_nodes,modifiedEdge,elemProperty] ...
                = bisection(splitEdgeFlags,...
                            node_coords,...
                            elem_nodes,...
                            edge_nodes,...
                            dualedge,...
                            elemProperty)   

    nNodes = size(node_coords,1);
    nElems = size(elem_nodes,1); 
    nEdges = size(edge_nodes,1);
    if (nargin < 5 || isempty(dualedge))
        dualedge = sparse(elem_nodes(:,[1,2,3]),...
                          elem_nodes(:,[2,3,1]),...
                          [1:nElems,1:nElems,1:nElems],...
                          nNodes,nNodes);
    end
    d2p = sparse(edge_nodes(:,[1,2]),...
                 edge_nodes(:,[2,1]),...
                 [1:nEdges,1:nEdges],...
                 nNodes,nNodes); 
    marker = zeros(nEdges,1);
    modifiedEdge= false(1,nEdges); % growing array is more efficient along last dimension

    nEdge = size(edge_nodes,1);
    %nEdgeUpdates = 0;
    newEdgeNum_nodes = [];
    % perform marking until a base is reached
    for base = 1:size(edge_nodes,1)
        if (splitEdgeFlags(base))
            flag = 1;
            el = dualedge(edge_nodes(base,1),edge_nodes(base,2));
            if (el == 0)
                el = dualedge(edge_nodes(base,2),edge_nodes(base,1));
            end
            while (flag)
                if (marker(base) > 0)
                    flag = 0;
                    % create new bisecting edge 
                    nEdge = nEdge+1;
                    newEdgeNum_nodes(:,end+1)  = [nEdge,elem_nodes(el,1),N]';
                    modifiedEdge(nEdge) = true;

                else
                    % mark midpoint of base and add the new node
                    N = size(node_coords,1)+1;
                    marker(base) = N;
                    node_coords(N,:) = mean(node_coords(elem_nodes(el,[2 3]),:));

                    % update edge and create new edge
                    newEdgeNum_nodes(:,end+1)  = [base,elem_nodes(el,2),N]'; % half of the original edge
                    modifiedEdge(base) = true;
                    nEdge = nEdge+1;
                    newSplitEdge = nEdge;
                    newEdgeNum_nodes(:,end+1)  = [nEdge,elem_nodes(el,3),N]'; % next half of the original edge
                    modifiedEdge(nEdge) = true;
                    nEdge = nEdge+1;
                    bisectEdge = nEdge;
                    newEdgeNum_nodes(:,end+1)  = [nEdge,elem_nodes(el,1),N]'; % new bisecting edge
                    modifiedEdge(nEdge) = true;
                    
                    right = d2p(elem_nodes(el,3),elem_nodes(el,1));
                    left = d2p(elem_nodes(el,1),elem_nodes(el,2));

                    if (marker(right) > 0)
                        % update edge and create new edge
                        nEdge = nEdge+1;
                        newEdgeNum_nodes(:,end+1) = [nEdge,sort([marker(base),marker(right)])]';
                        modifiedEdge(nEdge) = true;
                    end

                    if (marker(left) > 0)
                        % update edge and create new edge
                        nEdge = nEdge+1;
                        newEdgeNum_nodes(:,end+1) = [nEdge,sort([marker(base),marker(left)])]';
                        modifiedEdge(nEdge) = true;
                    end                    

                    % adjacent element that shares the same edge as a current
                    % element
                    el = dualedge(elem_nodes(el,3),elem_nodes(el,2));
                    
                    if (el == 0)
                        flag = 0;
                    else
                        elemMarkers(el) = false;
                    end
                end
            end
        end
    end
    
    if (nargin == 6)
        nProperty = length(elemProperty);
    end
    for i = 1:nElems
        base = d2p(elem_nodes(i,2),elem_nodes(i,3)); 
        if (marker(base) > 0)
           right = d2p(elem_nodes(i,3),elem_nodes(i,1));
           left = d2p(elem_nodes(i,1),elem_nodes(i,2));
           [elem_nodes,modifiedElems] ...
               = divide_elem(i,elem_nodes, ...
                             marker([base,right,left]));
           if (nargin == 6)
               for j = 1:nProperty
                    elemProperty{j}(modifiedElems) = elemProperty{j}(i);
               end
           end
        end
    end
    edge_nodes(newEdgeNum_nodes(1,:),:) = newEdgeNum_nodes(2:3,:)';

end



