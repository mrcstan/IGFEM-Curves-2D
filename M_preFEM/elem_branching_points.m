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
                                                               halfLineWidth,...
                                                               calcItrsectVel)
nElem = size(elem_nodes,1);
junc(1:nElem)=struct('nJunc',0);
nJunc = length(channels.junc);
newNode_coords = zeros(nJunc,2);
count = 0;
nEdgesPerElem = size(elem_edges,2);
for k = 1:nJunc
    xpt =channels.pts(channels.junc(k),1);
    ypt = channels.pts(channels.junc(k),2);
    for i = 1:nElem
        xElem = node_coords(elem_nodes(i,:),1);
        yElem = node_coords(elem_nodes(i,:),2);
        [IN,locEdge] = inTriangle(xpt,ypt,xElem,yElem,halfLineWidth);
        %IN =inpolygon(xpt,ypt,xElem,yElem);
        if (locEdge)  
            for j = 1:nEdgesPerElem
                [Lia,Locb] = ismemberf([xpt,ypt],...
                                       node_coords(itrsect(elem_edges(i,j)).node,:),...
                                       'row','tol',halfLineWidth);
                if (Lia)
                    break
                end
            end
			
			
            % if branching point lies on edge and has already been found by
            % edges_curves_intersect. In this case, the velocities
            % calculated by intersection_velocities.m is not correct and
            % has to be recalculated.
            if (Lia)
                 edge = elem_edges(i,locEdge);
                 if (calcItrsectVel)
                    itrsect(edge).vx{Locb} = [1,0];
                    itrsect(edge).vy{Locb} = [0,1];
                    [row,col] = find(channels.contvty==channels.junc(k),1,'first');
                    if(col(1) == 1)
                        ctrlPtNum = 1;
                    else
                        ctrlPtNum = channels.nurbs(row(1)).number;
                    end
                    itrsect(edge).designParamNum{Locb} ...
                        = channels.designParamNum{row(1)}(1:2,ctrlPtNum)';
                 end
            % if branching point lies on edge and has not be found by
            % edges_curves_intersect     
            else
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
                 if (calcItrsectVel)
                    itrsect(edge).vx{end+1} = [1,0];
                    itrsect(edge).vy{end+1} = [0,1];
                    
                    if(col(1) == 1)
                        ctrlPtNum = 1;
                    else
                        ctrlPtNum = channels.nurbs(row(1)).number;
                    end
                    itrsect(edge).designParamNum{end+1} ...
                        = channels.designParamNum{row(1)}(1:2,ctrlPtNum)';
                 end
            end
            break
        elseif (IN)
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
            break
        end
    end
end
%delrows = all(newNode_coords == 0,2);
%newNode_coords(delrows,:) = [];
newNode_coords = newNode_coords(1:count,:);
end



