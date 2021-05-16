%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Marcus Tan on 6/28/2014
%%% Last modified date: 6/28/2014
%%% Copyright 2014 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dualedge = elems_sharing_edge_nodes(elem_nodes,nNodes)
    repeatElems = repmat(1:size(elem_nodes,1),1,3);
    dualedge = sparse(elem_nodes(:,[1,2,3]),...
                      elem_nodes(:,[2,3,1]),...
                      repeatElems,...
                      nNodes,...
                      nNodes);
end