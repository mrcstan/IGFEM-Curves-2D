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
                = bisection(elemMarkers,...
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

    %nEdgeUpdates = 0;
    newEdgeNum_nodes = [];
    
    % the node of an element opposite its base
    elemPeaks = elem_nodes(:,1);
    
    elemDivided = false(nElems,1);
    % perform marking until a base is reached
    for i = 1:nElems
        el = i;
        if (elemMarkers(el))
            flag = 1;
            while (flag)
                base = d2p(elem_nodes(el,2),elem_nodes(el,3));
                if (marker(base) > 0)
                    flag = 0;                    
                    % create new bisecting edge 
                    nEdges = nEdges+1;
                    newEdgeNum_nodes(:,end+1)  = [nEdges,elemPeaks(el),N]';
                    modifiedEdge(nEdges) = true;
                    
                    
                    % need to update element peak because the peak of the element would
                    % have changed if it had been divided previously in this branch,  
                    % i.e, the element had been access via another element
                    % from the "else" block and now it is being access via
                    % yet another element from the "else" block
                    if (~elemDivided(el))
                        elemDivided(el) = true;
                        elemPeaks(el) = N;
                    end
                    
                    elemMarkers(el) = false;                   

                else
                    % mark midpoint of base and add the new node
                    N = size(node_coords,1)+1;
                    marker(base) = N;
                    node_coords(N,:) = mean(node_coords(elem_nodes(el,[2 3]),:));
                    
                   
                    
                    % update edge and create new edge
                    newEdgeNum_nodes(:,end+1)  = [base,elem_nodes(el,2),N]'; % half of the original edge
                    modifiedEdge(base) = true;
                    nEdges = nEdges+1;
                    newSplitEdge = nEdges;
                    newEdgeNum_nodes(:,end+1)  = [nEdges,elem_nodes(el,3),N]'; % next half of the original edge
                    modifiedEdge(nEdges) = true;
                    nEdges = nEdges+1;
                    bisectEdge = nEdges;
                    newEdgeNum_nodes(:,end+1)  = [nEdges,elem_nodes(el,1),N]'; % new bisecting edge
                    modifiedEdge(nEdges) = true;
                    
                    right = d2p(elem_nodes(el,3),elem_nodes(el,1));
                    left = d2p(elem_nodes(el,1),elem_nodes(el,2));

                    if (marker(right) > 0)
                        % update edge and create new edge
                        nEdges = nEdges+1;
                        newEdgeNum_nodes(:,end+1) = [nEdges,sort([marker(base),marker(right)])]';
                        modifiedEdge(nEdges) = true;
                    end

                    if (marker(left) > 0)
                        % update edge and create new edge
                        nEdges = nEdges+1;
                        newEdgeNum_nodes(:,end+1) = [nEdges,sort([marker(base),marker(left)])]';
                        modifiedEdge(nEdges) = true;
                    end                    
                    
                    elemDivided(el) = true;
                    % update element peak
                    elemPeaks(el) = N;
                    elemMarkers(el) = false;
                    
                    
                    % adjacent element that shares the same edge as a current
                    % element
                    el = dualedge(elem_nodes(el,3),elem_nodes(el,2));
                    
                    
                    if (el == 0)
                        flag = 0;

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



