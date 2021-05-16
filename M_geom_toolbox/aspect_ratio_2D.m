%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Marcus Tan on 10/4/2013
%%% Copyright 2013 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this function aspect ratio of a triangle or quadrilateral
% in 2D. It is defined as the max edge length/min edge length
%   INPUT:
%       coord: a 2x3 array [x1,x2,x3;y1,y2,y3] or [x1,x2,x3,x4;y1,y2,y3,y4];
function ratio = aspect_ratio_2D(coord)
edgeVec = diff(coord,1,2);
edgeVec(:,end+1) = coord(:,1)-coord(:,end);
edgeLens = sqrt(sum(edgeVec.^2,1));
ratio = max(edgeLens)/min(edgeLens);
end