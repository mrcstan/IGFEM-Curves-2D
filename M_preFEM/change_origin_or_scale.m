%%% Created by Marcus Tan on 1/31/2013
%%% Last modified date: 3/16/2014
%%% Copyright 2013 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this function translates the coordinates such that the new coordinates
% has their origin at Xo
% INPUT: 
%   nodeCoords
%   Xo: a length nDim array. nDim is the number of dimensions
%   Xnew = scale*(nodeCoords+Xo)
function nodeCoords = change_origin_or_scale(nodeCoords,Xo,scale)
nDim = size(nodeCoords,2);
for i = 1:nDim    
    nodeCoords(:,i) = scale*(nodeCoords(:,i)+Xo(i));
end
end