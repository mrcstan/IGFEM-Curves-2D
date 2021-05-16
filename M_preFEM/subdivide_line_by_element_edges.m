%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Marcus Tan on 1/22/2013
%%% Last modified date: 1/22/2014
%%% Copyright 2013 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subdivide the line using the tetrahedral element faces and further
% subdivide each segment into npts number of interior points
% INPUT:
%  line = [xo,yo,vx,vy] where (xo,yo) is the starting point
%  and (vx,vy) is the unit vector parallel to the line
%  npts: number of points between intersection with element edges
%  addLoc: additional points on the line described by the distance from 
%  (xo,yo)
function [xI,yI,sI,bm] = subdivide_line_by_element_edges(line,npts,nodeCoords,edge,addLoc)
tol0 = 1e-4;
tol1 = 1+tol0;
nEdges = size(edge.edge_node,1);
itrsectloc = [];
linePts = [line(1:2),line(1:2)+line(3:4)];
for i = 1:nEdges
    edgeLine = reshape(nodeCoords(edge.edge_node(i,:),:)',1,4);
    out = lineSegmentIntersect(edgeLine,linePts);
          
    if(out.intNormalizedDistance1To2 >-tol0 && out.intNormalizedDistance1To2<tol1)
        [~,~,itrsectloc] = find_and_update_sorted_real(out.intNormalizedDistance2To1,...
                                            itrsectloc,tol0);
    elseif(out.coincAdjacencyMatrix && edge.itrsect(i).num>0)
        for j = 1:edge.itrsect(i).num
            loc = linePosition3d([edge.itrsect(i).x(j),edge.itrsect(i).y(j),0],...
                                 [line(1:2),0,line(3:4),0]);
            [~,~,itrsectloc] = find_and_update_sorted_real(loc,itrsectloc,tol0);                 
        end
    end
end
if(nargin>4)
    for i = 1:length(addLoc)
        [~,~,itrsectloc] = find_and_update_sorted_real(addLoc(i),itrsectloc,tol0); 
    end
end
itrsectloc = sort(itrsectloc);
nItrsect = length(itrsectloc);
totNpts = (nItrsect-1)*npts+nItrsect;
xI = zeros(totNpts,1);
yI = zeros(totNpts,1);
sI = zeros(totNpts,1);
bm = false(totNpts,1);
for i=1:nItrsect-1
    sInd = (i-1)*(npts+1)+1;
    eInd = i*(npts+1)+1;
    s = linspace(itrsectloc(i),itrsectloc(i+1),npts+2);
    sI(sInd:eInd) = s;
    xI(sInd:eInd) = line(1)+s*line(3);
    yI(sInd:eInd) = line(2)+s*line(4);
    bm(sInd) = true;
    bm(eInd) = true;
end
end