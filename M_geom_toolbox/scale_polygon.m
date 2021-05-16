%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Marcus Tan on 4/3/2014
%%% Last modified date: 4/3/2014
%%% Copyright 2014 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this function expand or shrink a polygon
% verticesIn = [x1,x2,x3,...;y1,y2,y3,...]
function vertices = scale_polygon(verticesIn,scale)
vertices = verticesIn;
XXm = mean(vertices,2);
vertices(1,:) = vertices(1,:)-XXm(1);
vertices(2,:) = vertices(2,:)-XXm(2);
vertices = vertices*scale;
vertices(1,:) = vertices(1,:)+XXm(1);
vertices(2,:) = vertices(2,:)+XXm(2);
end