%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Marcus Tan on 4/21/2015
%%% Copyright 2015 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this function calculates the sparse matrix connecting the channel number
% and the vertex number to the control point number of the channel
% if i and j are the channel number and vertex number, respectively then
% the entry S(i,j) is the control point number
% INPUT:
%   option: 'c2v' = channel to vertex matrix
%           'v2c' = vertex to channel matrix
function [matrix,edgeList,nOrgVertices,nVertices] ...
                = channel2vertex_vertex2channel_matrix(channels,option)
nChannels = size(channels.contvty,1);
nT = sum(cat(1,channels.nurbs.number));
nIntVertices = nT - 2*nChannels;
nSmallEdges = nChannels + nIntVertices;
edgeList = nan(nSmallEdges,2);
edgeList(1:nChannels,:) = channels.contvty;
nOrgVertices = max(channels.contvty(:));
curVertex = nOrgVertices;
nVertices = curVertex + nIntVertices;
curEdge = nChannels;

ii = nan(nT,1);
jj = nan(nT,1);
ctrlPtNums = nan(nT,1);
ind = 1;
for i = 1:numel(channels.nurbs)
    ii(ind) = i;
    jj(ind) = channels.contvty(i,1);
    ctrlPtNums(ind) = 1; 
    ind = ind + 1;
    ii(ind) = i;
    jj(ind) = channels.contvty(i,2);
    ctrlPtNums(ind) = channels.nurbs(i).number;
    ind = ind + 1;
    if (channels.nurbs(i).number == 3)
        curEdge = curEdge + 1;
        edgeList(curEdge,2) = edgeList(i,2);
        curVertex = curVertex + 1;
        edgeList(i,2) = curVertex;
        edgeList(curEdge,1) = curVertex;
        ii(ind) = i;
        jj(ind) = curVertex;
        ctrlPtNums(ind) = 2;
        ind = ind + 1;
    elseif (channels.nurbs(i).number > 3)
        error('channels with more than 3 control points not yet considered')
    end
end
if (strcmpi(option,'c2v'))
    matrix = sparse(ii,jj,ctrlPtNums,nChannels,nVertices,nT);
elseif (strcmpi(option,'v2c'))
    matrix = sparse(jj,ii,ctrlPtNums,nVertices,nChannels,nT);
else
    error('unknown option')
end
end