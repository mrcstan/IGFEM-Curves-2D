%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Marcus Tan on 1/8/2015
%%% Copyright 2015 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [onEdges,normals] ...
                = kinks_branch_edges_intersect(KBpts, ...
                                               edges, ...
                                               tol)
onEdges = isPointOnEdge(KBpts,edges,tol);
if (size(KBpts,1) == 1)
    onEdges = find(onEdges);
else
    onEdges = find(any(onEdges,1));
end

if isempty(onEdges)
    normals = [];
    return
end

normals = [edges(onEdges,2) - edges(onEdges,4),...
           edges(onEdges,3) - edges(onEdges,1)];


normals = bsxfun(@rdivide,normals,sqrt(sum(normals.^2,2)));
end