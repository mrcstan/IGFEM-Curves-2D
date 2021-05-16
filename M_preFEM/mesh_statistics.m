%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Marcus Tan on 5/3/2014
%%% Copyright 2014 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stats = mesh_statistics(elem_nodes,nodeCoords)
maxval = -inf;
for i = 1:size(elem_nodes,1)
    cosMinAngle = cosine_min_angle_2D_triangle(nodeCoords(elem_nodes(i,:),:)');
    if (cosMinAngle > maxval)
        maxval = cosMinAngle;
    end
end
stats.minAngle = acos(maxval)*180/pi;
end