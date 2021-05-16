%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Marcus Tan on 7/24/2013
%%% Last modified date: 6/28/2014
%%% Copyright 2013 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this function determine the branching point% that are located in the
% interior and the edges (if it has not be found by edges_curves_intersect)
% of each element
% OUTPUT:
% junc.nJunc
% junc.pt
% node.n_node
% node.node_n
% node.coords
function [junc,newNode_coords,itrsect] = elem_branching_points(elem_nodes,...
                                                                 elem_edges,...
                                                                 itrsect,...
                                                                 node_coords,...
                                                                 nNodes,...
                                                                 channels,...
                                                                 halfLineWidth)
nElem = size(elem_nodes,1);
junc(1:nElem)=struct('nJunc',0);
nJunc = length(channels.junc);
newNode_coords = zeros(nJunc,2);
count = 0;
for k = 1:nJunc
    xpt =channels.pts(channels.junc(k),1);
    ypt = channels.pts(channels.junc(k),2);
    for i = 1:nElem
        xElem = node_coords(elem_nodes(i,:),1);
        yElem = node_coords(elem_nodes(i,:),2);
        [IN,locEdge] = inTriangle(xpt,ypt,xElem,yElem,halfLineWidth);
        %IN =inpolygon(xpt,ypt,xElem,yElem);
        if(IN && ~ismemberf([xpt,ypt],...
                           node_coords(cat(2,itrsect(elem_edges(i,:)).node),:),...
                           'row','tol',halfLineWidth)) 
            if (locEdge)   
                edge = elem_edges(i,locEdge);
                itrsect(edge).num = itrsect(edge).num+1;
                count = count+1;
                newNode_coords(count,:) = [xpt,ypt];  
                nNodes = nNodes+1;
                itrsect(edge).node(end+1) = nNodes;
                [row,col] = find(channels.contvty==channels.junc(k));
                itrsect(edge).nSeg(end+1) = length(row);
                itrsect(edge).seg{end+1} = row';
                itrsect(edge).x(end+1) = xpt;
                itrsect(edge).y(end+1) = ypt;
                itrsect(edge).nurbsParam{end+1} = col'-1;
            else                       
                junc(i).nJunc = junc(i).nJunc+1;
                junc(i).ipt(junc(i).nJunc) = channels.junc(k);
                [row,col] = find(channels.contvty==channels.junc(k));
                junc(i).nSeg(junc(i).nJunc) = length(row);
                junc(i).seg{junc(i).nJunc} = row';
                junc(i).nurbsParam{junc(i).nJunc} =col'-1; 
                count = count+1;  
                newNode_coords(count,:) = [xpt,ypt];  
                nNodes = nNodes+1;
                junc(i).node(junc(i).nJunc) = nNodes;
            end
            break;
        end
    end
end
delrows = all(newNode_coords == 0,2);
newNode_coords(delrows,:) = [];

end



