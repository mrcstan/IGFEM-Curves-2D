%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Marcus Tan on 12/22/2014
%%% Copyright 2014 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this function determines if a set of points fall on a set of NURBS curves
% within a given tolerance
% OUTPUT:
%   nodesOnNurbs: a vector of length nn containing the 
%                 global number of nodes coinciding with NURBS
%   normals: a nnx2 matrix, each row corresponding to the normal of the
%            NURBS on which the node lies on. if the node lies on the kink
%            of a NURBS consisting of straight line segments
%            the average normal direction of the two line segments meeting
%            at the kink is given
% REMARK: This is currently implemented for NURBS consisting of straight
%          segments only
function [nodesOnNurbs,normals] = nodes_curves_intersect(nodeCoords,nurbs,tolNode)
if (any(cat(1,nurbs.order) ~= 2))
    error('nodes_curves_intersect has not been implemented for NURBS curves')
end
nNurbs = numel(nurbs);
nEdges = sum(cat(1,nurbs.number)) - nNurbs;
edges = zeros(nEdges,4);
edgeNormals = zeros(nEdges,2);
count = 0;
for i = 1:nNurbs
    for j = 1:nurbs(i).number-1
        count = count + 1;
        edges(count,1:2) = nurbs(i).coefs(1:2,j)';
        edges(count,3:4) = nurbs(i).coefs(1:2,j+1)';
        edgeNormals(count,:) = [edges(count,2) - edges(count,4), ...
                                edges(count,3) - edges(count,1)];
        edgeNormals(count,:) = edgeNormals(count,:)/norm(edgeNormals(count,:));                    
    end
end
onNurbs = isPointOnEdge(nodeCoords,edges,tolNode);
if (size(nodeCoords,1) == 1)
    onNurbs = onNurbs';
end
nodesOnNurbs = find(any(onNurbs,2));

nn = numel(nodesOnNurbs);
normals = zeros(nn,2);
for i = 1:nn
    ind = onNurbs(nodesOnNurbs(i),:);
    normals(i,:) = mean(edgeNormals(ind,:),1);
    if (nnz(ind) > 1)
        normals(i,:) = normals(i,:)/norm(normals(i,:)); 
    end
end

end