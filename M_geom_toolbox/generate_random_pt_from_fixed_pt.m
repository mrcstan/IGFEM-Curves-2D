%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Marcus Tan on 3/21/2015
%%% Copyright 2015 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generage a random point a distance between lengthBounds(1) and
% lengthBounds(2) away from fixedPt and a polar angle between
% angleBounds(1) and angleBounds(2)
function [randPt,angle,dist] = generate_random_pt_from_fixed_pt(fixedPt, ...
                                                                lengthBounds, ...
                                                                angleBounds)
if (numel(fixedPt) ~= 2 || numel(lengthBounds) ~= 2 || numel(angleBounds) ~= 2)
    error('each input must be a length 2 vector')
end
angle = (angleBounds(2) - angleBounds(1))*rand() + angleBounds(1);
dist = (lengthBounds(2) - lengthBounds(1))*rand() + lengthBounds(1);
randPt(1) = fixedPt(1) + dist*cos(angle);
randPt(2) = fixedPt(2) + dist*sin(angle);
end