%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Marcus Tan on 11/12/2013
%%% Last modified date: 5/8/2014
%%% Copyright 2014 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This function finds the intersection of the edges of the elements with
% INPUT:
% edge_node
% node_coord
% elem_edge
% elem_region
% seg_region_region: a number of channels x 2 array. the i th row is a pair
%                   of region numbers between which the i th channel lies 
% OUTPUT:
% itrsect.num: number of intersections
% itrsect.nSeg: a length itrsect.num array specifying the number of line
%               source segments passing through a node 
% itrsect.seg: a itrsect.num x itrsect.nSeg array specifying the line 
%               source segments passing through a node 
% itrsect.pt: a length itrsect.num array indicating the terminal
%                     point of the line segment. if the intersection point
%                     is not a terminal point, it is set to zero.
% itrsect.node: a length itrsect.num array specifying the node number
% itrsect.x
% itrsect.y

function itrsect=edge_intersect_conforming(edge_node,node_coord, node_edges,...
                                            elem_edge,elem_region,itrface)
                                       
nEdge=size(edge_node,1);
itrsect = struct('num',0,'nSeg',cell(1,nEdge),'seg',cell(1,nEdge),...
                  'ipt',cell(1,nEdge),...
                  'node',cell(1,nEdge),'x',cell(1,nEdge),'y',cell(1,nEdge),...
                  'dist',cell(1,nEdge),'multItrsect',false);
itrsect = itrsect';
for i=1:nEdge % for i
    INd1=edge_node(i,1);
    INd2=edge_node(i,2);
    x1=node_coord(INd1,1);
    x2=node_coord(INd2,1);
    y1=node_coord(INd1,2);
    y2=node_coord(INd2,2);
    % find the elements sharing the edge given by ii
    [adjElem,~] = find(elem_edge == i);
    
    if(length(adjElem)==2)
        region = sort(elem_region(adjElem));
        if(region(1) ~= region(2))
           
            [row1,~] = find(itrface.seg_region_region == region(1));
            [row2,~] = find(itrface.seg_region_region == region(2));
            seg = intersect(row1,row2,'stable');
            % decides which channel is the right one if multiple channels
            % are in between pairs of the same region
            if (length(seg)>1)
                for j = 1:length(seg)
                    for k = 1:itrface.nurbs(seg(j)).nConvHull
                        if(itrface.nurbs(seg(j)).iscollinear{k})
                            xmin = min(itrface.nurbs(seg(j)).hull{k}(1),itrface.nurbs(seg(j)).hull{k}(3));
                            ymin = min(itrface.nurbs(seg(j)).hull{k}(2),itrface.nurbs(seg(j)).hull{k}(4));
                            xmax = max(itrface.nurbs(seg(j)).hull{k}(1),itrface.nurbs(seg(j)).hull{k}(3));
                            ymax = max(itrface.nurbs(seg(j)).hull{k}(2),itrface.nurbs(seg(j)).hull{k}(4));
                            hull = [xmin,xmax,xmax,xmin;ymin,ymin,ymax,ymax];
                           
                        else
                            hull = itrface.nurbs(seg(j)).hull{k};
                        end
                        intersectHull = isintersect(hull,[x1,x2;y1,y2]);
                        %[~,intersectHull] = intersectEdgePolygon([x1,y1,x2,y2], hull');
                        if (intersectHull)
                            break
                        end
                    end
                    if (intersectHull)
                        seg = seg(j);
                        break
                        
                    end
                end
                if (length(seg)>1)
                    error('none of the convex hulls contain the channels')
                end
            end
            itrsect(i).num = 2;
            if(iscell(itrsect(i).seg))
                itrsect(i).seg{1}(end+1) = seg;                
                itrsect(i).seg{2}(end+1) = seg;
                itrsect(i).nSeg(1) = itrsect(i).nSeg(1)+1;
                itrsect(i).nSeg(2) = itrsect(i).nSeg(2)+1;
            else
                itrsect(i).seg{1} = seg;
                itrsect(i).seg{2} = seg;
                itrsect(i).nSeg(1) = 1;
                itrsect(i).nSeg(2) = 1;
            end            
            itrsect(i).pt(1:2) = 0;
            itrsect(i).node(1:2) = [INd1,INd2];
            itrsect(i).x(1:2) = [x1,x2];
            itrsect(i).y(1:2) = [y1,y2];
            
            % edge node 1
            sharedEdges = node_edges{INd1};
            for j = 1:length(sharedEdges)
                edge = sharedEdges(j);
                ind = find(itrsect(edge).node == INd1,1,'first');
                if (~isempty(ind) && all(itrsect(edge).seg{ind} ~= seg))
                    itrsect(edge).nSeg(ind) = itrsect(edge).nSeg(ind)+1;
                    itrsect(edge).seg{ind} = [itrsect(edge).seg{ind},seg];
                end
            end
            
            % edge node2
            sharedEdges = node_edges{INd2};
            for j = 1:length(sharedEdges)
                edge = sharedEdges(j);
                ind = find(itrsect(edge).node == INd2,1,'first');
                if (~isempty(ind) && all(itrsect(edge).seg{ind} ~= seg))
                    itrsect(edge).nSeg(ind) = itrsect(edge).nSeg(ind)+1;
                    itrsect(edge).seg{ind} = [itrsect(edge).seg{ind},seg];
                end
            end
        end
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Marcus Tan on 7/23/2013
%%% Last modified date: 7/23/2013
%%% Copyright 2013 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% this function finds all segments sharing a point
function [seg,nSeg,pt]=find_seg_sharing_pt(ipt,contvty)
[row,col]=find(contvty==ipt);
nSeg=length(row);
seg=row';
if(nSeg>1)
    pt=ipt;
else
    pt=0;
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Marcus Tan on 7/24/2013
%%% Last modified date: 7/24/2013
%%% Copyright 2013 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% this function determines whether an intersection point which is a 
% terminal point of a line segment has already been found.
function flag=pt_considered(pt,itrsect)
flag=false;
if(isfield(itrsect,'pt'))
    if(any(itrsect.pt==pt))
        flag=true;
    end
end
end





