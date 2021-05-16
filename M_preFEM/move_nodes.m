%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Marcus Tan on 12/22/2014
%%% Modified on 1/16/2016
%%% Copyright 2014 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this function moves a given set of nodes by a distance del either in
% a fixed direction or randomly
% INPUT:
%   nodeCoords: nn x 2 matrix, each row contains the coordinates of each
%               node to be moved 
%   directions: nn x 2 matrix, each row contains the direction each node
%               should move
function nodeCoords = move_nodes(nodeCoords,del,directions,randDirFlag)

if randDirFlag
    dirSign = sign(rand(size(nodeCoords,1),1)-0.5);
else
    dirSign = 1.0;
end
nodeCoords = nodeCoords + del*bsxfun(@times,directions,dirSign);


end