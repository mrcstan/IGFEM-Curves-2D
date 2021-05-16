%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Marcus Tan on 3/21/2015
%%% Copyright 2015 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this function generates a random quadrilateral
% INPUT:
%   coordBounds: [xmin,xmax;ymin,ymax] = bounds for vertex coordinates 
%   lengthBounds: [r1min,r1max;r2min,r2max;r3min,r3max;r4min,r4max] 
%                   = bounds on the lengths of edges 12, 23, 34 and 14
%                some of the bounds may be immaterial if the vertices are
%                already fixed
%   minAngle: min interior angle. max allowed interior angle would be pi -
%             minAngle
%   fixedPoints: [x1,y1;...;xm,ym] = fixed vertices in CCW order, 1<=m<4
%                specifiy [] if no fixed point
%   maxTrials: maximum number of attempts
%   showWarning: show warning when generation is unsuccessful
% OUTPUT:
%   quad: [x1,y1;x2,y2;x3,y3;x4,y4]
%   generated: true if quad satisfying specifications is successfully
%               generated
%   attempts: actual number of attempts
function [quad,generated,attempts] = generate_random_quad(coordBounds, ...
                                                lengthBounds, ...
                                                minAngle, ...
                                                fixedPoints, ...
                                                maxTrials,...
                                                showWarning)
if (size(coordBounds,1) ~= 2 || size(coordBounds,2) ~= 2)
    error('coordBounds should be a 2x2 matrix')
end
if (size(lengthBounds,2) ~= 2)
    error('lengthBounds must have 2 columns')
end

if (size(lengthBounds,1) ~= 4)
    error('lengthBounds must have 4 rows')
end
if (minAngle <= 0 || minAngle >= 0.5*pi)
    error('minAngle must be between 0 and pi/2')
else
    sinMinAngle = sin(minAngle);
    maxAngle = pi - minAngle;
end
nFixedPts = size(fixedPoints,1);
if (nFixedPts > 3)
    error('number of fixed points < 4')
end

quad = nan(4,2);
if (nFixedPts)
    quad(1:nFixedPts,:) = fixedPoints;
end
r4minsq = lengthBounds(4,1)^2;
r4maxsq = lengthBounds(4,2)^2;
delx = coordBounds(1,2) - coordBounds(1,1);
dely = coordBounds(2,2) - coordBounds(2,1);
delLengths = lengthBounds(:,2) - lengthBounds(:,1);
delAngle = maxAngle - minAngle;
generated = false;
rng('shuffle')

if (nFixedPts == 2)
    [angle1,dist1] = cart2pol(quad(2,1)-quad(1,1),quad(2,2)-quad(1,2));
    if (dist1 < lengthBounds(1,1) || dist1 > lengthBounds(1,2))
        error('the distance between first two vertices is out of bound')
    end
end
%
if (nFixedPts == 3)
    [angle3,dist3] = cart2pol(quad(3,1)-quad(2,1),quad(3,2)-quad(2,2));
    if (dist3 < lengthBounds(2,1) || dist3 > lengthBounds(2,2))
        error('the distance between 2nd and 3rd vertices is out of bound')
    end
end
% generate first point
if (nFixedPts == 0)
    quad(1,1) = delx*rand() - coordBounds(1,1);
    quad(1,2) = dely*rand() - coordBounds(2,1);
end
    
for attempts = 1:maxTrials
   
    % generate 2nd point
    if (nFixedPts <= 1)
        %maxAngle1 = 0.5*pi - minAngle;
        [quad(2,:),angle1] = generate_random_pt_from_fixed_pt(quad(1,:), ...
                                                              lengthBounds(1,:), ...
                                                              [-pi,pi]);
    end
    % check if 2nd point is out of bound
    if (~in_box(quad(2,:),coordBounds))
        continue
    end
    
    if (nFixedPts <= 2)
        angle2 = delAngle*rand() + minAngle;
        randDist = delLengths(2)*rand() + lengthBounds(2,1);
        angle3 = angle1 + angle2;
        quad(3,1) = quad(2,1) + randDist*cos(angle3);
        quad(3,2) = quad(2,2) + randDist*sin(angle3);
    end
    % check if 3rd point is out of bound
    if (~in_box(quad(3,:),coordBounds))
        continue
    end
    %
    if (nFixedPts <= 3)
        angle4 = delAngle*rand() + minAngle;
        randDist = delLengths(3)*rand() + lengthBounds(3,1);
        angle5 = angle3 + angle4;
        quad(4,1) = quad(3,1) + randDist*cos(angle5);
        quad(4,2) = quad(3,2) + randDist*sin(angle5);
    end
    % check if 4th point is out of bound
    if (~in_box(quad(4,:),coordBounds))
        continue
    end
    
    dist14sq = (quad(4,1)-quad(1,1))^2 + (quad(4,2)-quad(1,2))^2;
    if (dist14sq < r4minsq || dist14sq > r4maxsq)
        continue
    end
    
    sinAngles = polygon_sine_angles(quad(:,1),quad(:,2));
    if (any(sinAngles < sinMinAngle))
        continue
    end
    
    
    generated = true;
    break
    %{
    dist13sq = (quad(3,1)-quad(1,1))^2 + (quad(3,2)-quad(1,2))^2;
    if (dist13sq > 0.25*rmaxsq)
        continue;
    end
    if (nFixedPts <= 3)
        minDist = 0.5*sqrt(dist13sq);
        randDist = (lengthBounds(2) - minDist)*rand(1,2) + minDist;
        [xout,yout] = circcirc(quad(1,1),quad(1,2),randDist(1),...
                                         quad(3,1),quad(3,2),randDist(2));       
    end
    if (isnan(xout))
        continue
    end
    for i = 1:numel(xout)
        quad(4,1) = xout(i);
        quad(4,2) = yout(i);
        sinAngles = polygon_sine_angles(quad(:,1),quad(:,2));
        if (all(sinAngles > sinMinAngle))
            generated = true;
            break
        end
    end
    if(generated)
        break
    end
    %}
    
end

if (~generated && showWarning)
    warning('maximum attempts exceed, generated quad does not satisfy specs')
    %quad = [];
end
end