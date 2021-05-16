%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Marcus Tan on 10/4/2013
%%% Copyright 2013 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this function calculates the cosine of the smallest angle in a triangle
% in 2D
%   INPUT:
%       coord: a 2x3 array [x1,x2,x3;y1,y2,y3];
function cosMinAngle = cosine_min_angle_2D_triangle(coord)
lensq = zeros(3,1);
ind = [1,2;
       2,3;
       3,1];
for i=1:3
    vec = coord(:,ind(i,2))-coord(:,ind(i,1));
    lensq(i) = dot(vec,vec);
end
cosAngle = -10.0*ones(3,1);
ind = [2,3,1;
       1,3,2;
       1,2,3];
for i=1:3
    cosAngle(i) = 0.5*(lensq(ind(i,1))+lensq(ind(i,2))-lensq(ind(i,3)))/(sqrt(lensq(ind(i,1))*lensq(ind(i,2)))); 
end
cosMinAngle = max(cosAngle);
end