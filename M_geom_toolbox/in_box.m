%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Marcus Tan on 3/21/2015
%%% Copyright 2015 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this function determines if a point is in a box
% INPUT:
%   boxBounds: [xmin,xmax;ymin,ymax]
%   pt: [x,y]
function in = in_box(pt,boxBounds)
if (pt(1) < boxBounds(1,1) || pt(1) > boxBounds(1,2) ...
  ||pt(2) < boxBounds(2,1) || pt(2) > boxBounds(2,2))
    in = false;
else
    in = true;
end
end