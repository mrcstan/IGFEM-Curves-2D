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
                                                        halfLineWidth, ...
                                                        calcItrsectVel, ...
                                                        channelDesignParamNum)
nElem = size(elem_nodes,1);
elemKinks(1:nElem)=struct('nKinks',0);
nKinks = length(kinks.x);
newNode_coords = zeros(nKinks,2);
count = 0;
nEdgesPerElem = size(elem_edges,2);
for k = 1:nKinks
    xpt = kinks.x(k);
    ypt = kinks.y(k);
    for i = 1:nElem
        xElem = node_coords(elem_nodes(i,:),1);
        yElem = node_coords(elem_nodes(i,:),2);
        [IN,locEdge] = inTriangle(xpt,ypt,xElem,yElem,halfLineWidth);

        %[IN,ON] = inpolygon(xpt,ypt,xElem,yElem);
        if(locEdge)
            for j = 1:nEdgesPerElem
                [Lia,Locb] = ismemberf([xpt,ypt],...
                                       node_coords(itrsect(elem_edges(i,j)).node,:),...
                                       'row','tol',halfLineWidth);
                if (Lia)
                    break
                end
            end
            % if kink lies on edge and has already been found by
            % edges_curves_intersect. In this case, the velocities
            % calculated by intersection_velocities.m is not correct and
            % has to be recalculated.
            if (Lia)
                if (calcItrsectVel)
                    edge = elem_edges(i,locEdge);
                    itrsect(edge).vx{Locb} = [1,0];
                    itrsect(edge).vy{Locb} = [0,1];
                    itrsect(edge).designParamNum{Locb} ...
                        = channelDesignParamNum{kinks.seg(k)}(1:2,kinks.ctrlPtNum(k))';
                end
            % if kink lies on edge and has not be found by
            % edges_curves_intersect    
            else
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
                if (calcItrsectVel)
                    itrsect(edge).vx{end+1} = [1,0];
                    itrsect(edge).vy{end+1} = [0,1];
                    itrsect(edge).designParamNum{end+1} ...
                        = channelDesignParamNum{kinks.seg(k)}(1:2,kinks.ctrlPtNum(k))';
                end
            end
            break
        elseif(IN)
            elemKinks(i).nKinks = elemKinks(i).nKinks+1;
            elemKinks(i).nSeg(elemKinks(i).nKinks) = 1;
            elemKinks(i).seg{elemKinks(i).nKinks} = kinks.seg(k);
            elemKinks(i).nurbsParam{elemKinks(i).nKinks} = kinks.nurbsParam(k);
            elemKinks(i).ctrlPtNum(elemKinks(i).nKinks) = kinks.ctrlPtNum(k);
            count = count+1;
            newNode_coords(count,:) = [xpt,ypt];  
            nNodes = nNodes+1;
            elemKinks(i).node(elemKinks(i).nKinks) = nNodes;
            break
        end
    end
end
%delrows = all(newNode_coords == 0,2);
%newNode_coords(delrows,:) = [];
newNode_coords = newNode_coords(1:count,:);

end

