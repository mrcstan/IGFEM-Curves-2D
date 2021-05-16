%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Marcus Tan on 11/12/2013
%%% Last modified date: 11/12/2013
%%% Copyright 2013 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This function finds the intersection of the edges of the elements with
% INPUT:
% edge_node
% node_coord
% elem_edge
% elem_region
% region_region_seg: a sparse matrix specifying the interface number
%                    between two regions
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

function itrsect=edge_intersect_conforming(edge_node,node_coord,...
                                             elem_edge,elem_region,region_region_seg)
nEdge=size(edge_node,1);
itrsect = struct('num',0,'nSeg',cell(1,nEdge),'seg',cell(1,nEdge),...
                  'ipt',cell(1,nEdge),...
                  'node',cell(1,nEdge),'x',cell(1,nEdge),'y',cell(1,nEdge),...
                  'dist',cell(1,nEdge),'doubleItrsect',false);
for i=1:nEdge % for i
    INd1=edge_node(i,1);
    INd2=edge_node(i,2);
    x1=node_coord(INd1,1);
    x2=node_coord(INd2,1);
    y1=node_coord(INd1,2);
    y2=node_coord(INd2,2);
    % find the elements sharing the edge given by ii
    adjElem1 = find(elem_edge(:,1)==i);
    adjElem2 = find(elem_edge(:,2)==i);
    adjElem3 = find(elem_edge(:,3)==i);
    adjElem = [adjElem1',adjElem2',adjElem3'];
    if(length(adjElem)==2)
        region = sort(elem_region(adjElem));
        if(region(1) ~= region(2))
            itrsect(i).num=2;
            itrsect(i).nSeg(1:2) = 1;
            itrsect(i).seg{1} = region_region_seg(region(1),region(2));
            if(isempty(itrsect(i).seg{1}))
                warning('cannot assign a label to interface between region') 
            end
            itrsect(i).seg{2} = itrsect(i).seg{1};
            itrsect(i).pt(1:2) = 0;
            itrsect(i).node(1:2) = [INd1,INd2];
            itrsect(i).x(1:2) = [x1,x2];
            itrsect(i).y(1:2) = [y1,y2];
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Marcus Tan on 7/25/2013
%%% Last modified date: 7/25/2013
%%% Copyright 2013 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function itrsectNew=update_itrsect(itrsectOld,seg,nSeg,pt,node,nodeCoords)
    itrsectNew=itrsectOld;
    itrsectNew.num=itrsectNew.num+1;
    num=itrsectNew.num;
    itrsectNew.nSeg(num)=nSeg;
    itrsectNew.seg{num}=seg;
    itrsectNew.pt(num)=pt;
    itrsectNew.node(num)=node;
    itrsectNew.x(num)=nodeCoords(1);
    itrsectNew.y(num)=nodeCoords(2);   
end




