%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Marcus Tan on 6/12/2013
%%% Last modified date: 3/15/2013
%%% Copyright 2013 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This function added interpolate temperature value to the nodal value of
% new added intersection nodes.

function UUR2 = update_enrichment_node_value(UUR,node,edge,elem)

UUR2=UUR;
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
    if(~isempty(elem.junc))
        if(elem.junc(i).nJunc>0)
            shape=1;
            for j=1:elem.junc(i).nJunc
                Xel=[node.coords(elem.elem_node(i,:),1),...
                     node.coords(elem.elem_node(i,:),2)]';
                curNode=elem.junc(i).node(j);
                Xlocal=local_coord(shape,node.coords(curNode,:),Xel);
               [N,~]=shape_funct(Xlocal, shape);    
               UUR2(curNode)=UUR(elem.elem_node(i,:))'*N+UUR(curNode);
            end
        end
    end
    if(elem.parent(i).isParent)
        shape=1;
        Xel=[node.coords(elem.elem_node(i,:),1),...
                 node.coords(elem.elem_node(i,:),2)]';
        for j=1:elem.parent(i).nNode
            curNode=elem.parent(i).node(j);
            if(curNode>node.nOriginalEnrichNode)
                Xlocal=local_coord(shape,node.coords(curNode,:),Xel);
                [N,~]=shape_funct(Xlocal, shape);   
                UUR2(curNode)=UUR(elem.elem_node(i,:))'*N+UUR(curNode);             
            end
        end
    end
end

