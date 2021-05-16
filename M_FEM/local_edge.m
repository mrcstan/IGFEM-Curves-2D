%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Marcus Tan on 7/16/2013
%%% Last modified date: 7/16/2013
%%% Copyright 2013 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% this function returns the local edge number of the current edge defined
% by the two nodes in the array inNode
% INPUT:
% refNode: a length 3 or 4 array with the nodes ordered counter clockwise
% inNode: a length 2 array defining the current edge
% OUTPUT:
% locEdge: local edge number
% outNode: a length 2 array such that the vector from outNode(1) to
% outNode(2) points in the direction of increasing natural coordinate
function [locEdge,outNode]=local_edge(refNode,inNode)
nRefNode=length(refNode);
outNode=zeros(2,1);
node1=inNode(1);
node2=inNode(2);
i1=find(refNode==node1);
i2=find(refNode==node2);

if(isempty(i1)||isempty(i2))
    error('local_edge: at least one of the nodes in inNode is not found in refNode');
end
if(i1==4 && i2==3)

elseif(i1>i2 || (i1==3 && i2==4))
    i3=i2;
    temp=node2;
    i2=i1;
    node2=node1;
    i1=i3;
    node1=temp;
end
if(nRefNode==4)
    if(i1==1 && i2==2)
        locEdge=1;
    elseif(i1==2 && i2==3)
        locEdge=2;
    elseif(i1==4 && i2==3)
        locEdge=3;
    elseif(i1==1 && i2==4)
        locEdge=4;
    end
elseif(nRefNode==3)
    if(i1==1 && i2==2)
        locEdge=1;
    elseif(i1==2 && i2==3)
        locEdge=2;
    elseif(i1==1 && i2==3)
        locEdge=3;
    end
else
    disp('local_edge: length of refNode should be 3 or 4')
end
outNode(1)=node1;
outNode(2)=node2;
end