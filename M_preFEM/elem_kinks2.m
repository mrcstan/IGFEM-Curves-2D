%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Marcus Tan on 7/24/2013
%%% Last modified date: 7/11/2014
%%% Copyright 2013 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this function determine the kinks of a channel
% that are located in the interior and the edges (if it has not be found by 
% edges_curves_intersect) of each element. Kinks are only present
% for a channel consisting of straight line segments
% OUTPUT:
% elemKinks(i).nKinks
% elemKinks(i).nSeg
% elemKinks(i).seg
% elemKinks(i).nurbsParam
% elemKinks(i).ctrlPtNum
% newNode_coords
% itrsect
function [elemKinks,newNode_coords,itrsect] =elem_kinks(elem_nodes, ...
                                                        elem_edges, ...
                                                        itrsect, ...
                                                        node_coords, ...
                                                        nNodes,kinks, ...
                                                        halfLineWidth)
nElem = size(elem_nodes,1);
elemKinks(1:nElem)=struct('nKinks',0);
nKinks = length(kinks.x);
newNode_coords = zeros(nKinks,2);
count = 0;
for k = 1:nKinks
    xpt = kinks.x(k);
    ypt = kinks.y(k);
    for i = 1:nElem
        xElem = node_coords(elem_nodes(i,:),1);
        yElem = node_coords(elem_nodes(i,:),2);
        [IN,locEdge] = inTriangle(xpt,ypt,xElem,yElem,halfLineWidth);

        %[IN,ON] = inpolygon(xpt,ypt,xElem,yElem);
        if(IN && ~ismemberf([xpt,ypt],...
                             node_coords(cat(2,itrsect(elem_edges(i,:)).node),:),...
                             'row','tol',halfLineWidth))
            % if kink lies on edge and has not be found by
            % edges_curves_intersect
            if (locEdge) 
                edge = elem_edges(i,locEdge);
                itrsect(edge).num = itrsect(edge).num+1;
                count = count+1;
                newNode_coords(count,:) = [xpt,ypt];  
                nNodes = nNodes+1;
                itrsect(edge).node(end+1) = nNodes;
                itrsect(edge).nSeg(end+1) = 1;
                itrsect(edge).seg{end+1} = kinks.seg(k);
                itrsect(edge).x(end+1) = xpt;
                itrsect(edge).y(end+1) = ypt;
                itrsect(edge).nurbsParam{end+1} = kinks.nurbsParam(k);
            else
                elemKinks(i).nKinks = elemKinks(i).nKinks+1;
                elemKinks(i).nSeg(elemKinks(i).nKinks) = 1;
                elemKinks(i).seg{elemKinks(i).nKinks} = kinks.seg(k);
                elemKinks(i).nurbsParam{elemKinks(i).nKinks} = kinks.nurbsParam(k);
                elemKinks(i).ctrlPtNum(elemKinks(i).nKinks) = kinks.ctrlPtNum(k);
                count = count+1;
                newNode_coords(count,:) = [xpt,ypt];  
                nNodes = nNodes+1;
                elemKinks(i).node(elemKinks(i).nKinks) = nNodes;
            end
            break
        end
    end
end
delrows = all(newNode_coords == 0,2);
newNode_coords(delrows,:) = [];

end

