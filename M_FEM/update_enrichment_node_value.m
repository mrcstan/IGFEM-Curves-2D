%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Marcus Tan on 6/12/2013
%%% Last modified date: 3/15/2013
%%% Copyright 2013 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This function added interpolate temperature value to the nodal value of
% new added intersection nodes.
% INPUT: 
%   epsge: geometric tolerance for nurbsSurf_phy2param. Default is 1e-8.
function UUR2 = update_enrichment_node_value(UUR,node,elem,edge,epsge)

UUR2 = UUR;
shape = 1;
if (nargin > 4)
    % evaluate the true enrichment node values (slow)
        ind = (node.nOriginalNode+1):node.n_node;
        [u_du,~,~,del] = interpolate_soln(node.coords, ...
                                            elem, ...
                                            UUR, ...
                                            node.coords(ind,1),...
                                            node.coords(ind,2),...
                                            epsge,...
                                            false);
        del = del(del > 0) + node.nOriginalNode;
        UUR2(setdiff(ind,del)) = u_du(:,1);
        if (isempty(del))
            return
        else
            fprintf('%i enrichment nodes cannot be evaluate \n',length(del))
        end

     % approximate the true enrichment node values with linear interpolations (fast)
else
        nEdge=size(edge.edge_node,1);
        % update temperature value at enrichment nodes on edge
        for i=1:nEdge
            if(edge.itrsect(i).num>0)
                for j=1:edge.itrsect(i).num
                    curNode=edge.itrsect(i).node(j);
                    if(curNode>node.nOriginalNode)
                         endNode1=edge.edge_node(i,1);
                         endNode2=edge.edge_node(i,2);
                         d=sqrt((node.coords(endNode2,1)-node.coords(endNode1,1))^2....
                               +(node.coords(endNode2,2)-node.coords(endNode1,2))^2);
                         d1=sqrt((node.coords(curNode,1)-node.coords(endNode1,1))^2....
                              +(node.coords(curNode,2)-node.coords(endNode1,2))^2);
                         UUR2(curNode) = UUR(curNode)+UUR(endNode1)...
                                        +d1/d*(UUR(endNode2)-UUR(endNode1));
                    end
                end
            end
        end
        % update temperature value at enrichment node in interior of element
        for i=1:elem.n_elem
            if ~isempty(elem.branch_kinks(i).bkNums)
                for j = 1:numel(elem.branch_kinks(i).bkNums)
                    Xel = node.coords(elem.elem_node(i,:),:)';
                    curNode = elem.branch_kinks(i).nodes(j);
                    Xlocal = local_coord(node.coords(curNode,:)',Xel,1); % assume parent element is always a triangular element
                   [N,~] = shape_funct(Xlocal, shape);    
                   UUR2(curNode) = UUR(elem.elem_node(i,:))'*N+UUR(curNode);
                end
            end            
            if(elem.parent(i).type > 1)
                Xel=[node.coords(elem.elem_node(i,:),1),...
                         node.coords(elem.elem_node(i,:),2)]';
                for j=1:numel(elem.parent(i).nodes)
                    curNode=elem.parent(i).nodes(j);
                    if(curNode>node.nOriginalEnrichNode)
                        Xlocal=local_coord(node.coords(curNode,:)',Xel,1); % assume parent element is always a triangular element
                        [N,~]=shape_funct(Xlocal, shape);   
                        UUR2(curNode)=UUR(elem.elem_node(i,:))'*N+UUR(curNode);             
                    end
                end
            end
        end
end
end

