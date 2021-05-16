%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Marcus Tan on 3/6/2016
%%% Copyright 2016 University of Illinois at Urbana-Champaign. 
%%% All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Determine the local edge of a triangle given the barycentric 
% coordinates of a point in the triangle. The barycentric coordinates
% are defined by MATLAB pointLocation function
% If the point is in the interior, the local face is defined as 0
% INPUT:
%   locCoords: a vector with 3 entries
% OUTPUT:
%   locFace: the local face numbering should be consistent with that in
%            local_face.m
function [locEdge,edgeLocNodes] ...
                = loc_edge_from_bary_coords(baryCoords,tol)
if numel(baryCoords) ~= 3
    error('barycentric coordinates must have 3 entries')
end
allEdgeLocNodes = [1,2;
                   2,3;
                   1,3]';
locEdge = 0;
edgeLocNodes = zeros(2,1);

for i = 1:size(allEdgeLocNodes,2);
    if (abs(sum(baryCoords(allEdgeLocNodes(:,i)))-1) < tol)
        locEdge = i;
        edgeLocNodes = allEdgeLocNodes(:,i);
    end
end
end