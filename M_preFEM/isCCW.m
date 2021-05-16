%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Marcus Tan on 8/25/13
%%% Modified by Marcus Tan
%%% Last modified date: 7/15/14
%%% Copyright 2013 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this function determines if a pair of local nodes are counter-clockwise
% or not. It is assumed that the counter-clockwise reference nodes are
% labelled [1,2,...,maxLocNode]
% INPUT: 
%   locNode: length 2 array
function flag=isCCW(locNode,maxLocNode)
flag=false;
if(locNode(1)==1 && locNode(2)==maxLocNode)
    flag=false;
elseif(locNode(1)==maxLocNode && locNode(2)==1)
    flag=true;
elseif(locNode(2)>locNode(1))
    flag=true;
end
end