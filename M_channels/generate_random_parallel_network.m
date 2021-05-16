%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Marcus Tan on 3/22/2015
%%% Copyright 2015 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function generates a parallel network based on the number of
% branches specified
% INPUT:
%   iniVertex: inlet vertex
%   outVertex: outlet vertex
%   nBranches: number of branches
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
%   vertices: a 2 column matrix of the coordinates of the vertices
%   connectivity: a (nBranches+2)x4 connectivity matrix of the
%                 quadrilaterals forming the loops of the network and 
%                 the two initial and final line segments
%                 only first two entries of the first and last rows are
%                 meaningful since a line segment connects only 2 vertices
%   generated: true if network successfully generated
function [vertices,connectivity,generated] ...
                = generate_random_parallel_network(iniVertex, ...
                                                   finalVertex, ...
                                                   nBranches, ...
                                                   coordBounds, ...
                                                   lengthBounds, ...
                                                   minAngle, ...
                                                   maxTrials, ...
                                                   showWarning)
if (nBranches < 2)
    error('number of branches must be greater than 2')
end
if (numel(iniVertex) ~= 2 || numel(finalVertex) ~= 2)
    error('iniVertex,finalVertex,lengthBounds must have 2 entries')
end
if (size(coordBounds,1) ~= 2 || size(coordBounds,2) ~= 2)
    error('coordBounds must be a 2x2 matrix')
end

if (size(lengthBounds,2) ~= 2)
    error('lengthBounds must have 2 columns')
end

if (size(lengthBounds,1) ~= 4)
    error('lengthBounds must have 4 rows')
end

if (maxTrials < 1)
    error('maxTrials >= 1')
end

firstLengthBounds = circshift(lengthBounds,[-1,0]); % for first loop
nVertices = 2*(nBranches + 2);
vertices = nan(nVertices,2);
vertices(1,:) = iniVertex;
vertices(end,:) = finalVertex;
maxRandQuadTrials = 200;
generated = false;

% create connectivity matrix
%nChannels = 3*nBranches + 1;

connectivity = nan(nBranches+2,4);
connectivity(1,1:2) = [1,2];
connectivity(end,1:2) = [nVertices - 1,nVertices];
for i = 1:nBranches
    iTopLeft = 2*i;
    connectivity(i+1,:) = [iTopLeft, iTopLeft+2, iTopLeft+3, iTopLeft+1];
end


for attempts = 1:maxTrials
    vertices(2,:) = generate_random_pt_from_fixed_pt(iniVertex, ...
                                                     [lengthBounds(2,1),lengthBounds(1,2)], ...
                                                     [-pi/2,pi/2]);
   if (~in_box(vertices(2,:),coordBounds))
       continue
   end
   
   % create first loop
   [quad,quadGenerated] = generate_random_quad(coordBounds, ...
                                               firstLengthBounds, ...
                                               minAngle, ...
                                               vertices(2,:), ...
                                               maxRandQuadTrials, ...
                                               false);
   if (~quadGenerated)
       continue
   end
   
   % check if first line segment intersects with the first quad
   xi = polyxpoly(vertices(1:2,1),vertices(1:2,2),quad(:,1),quad(:,2));
   if (numel(xi) > 1)
       %disp('first segment intersects first quad')
       continue
   end
   vertices(3,:) = quad(4,:);
   vertices(4,:) = quad(2,:);
   vertices(5,:) = quad(3,:);
   quad = vertices(2:5,:);   
   % create remaining loops
   for i = 1:(nBranches - 1)
        [quad,quadGenerated] = generate_random_quad(coordBounds, ...
                                                   lengthBounds, ...
                                                   minAngle, ...
                                                   quad(4:-1:3,:), ...
                                                   maxRandQuadTrials, ...
                                                   false);
                                               
        if (~quadGenerated)
            break
        end            
        % check if first line segment intersects with the first quad
        xi = polyxpoly(vertices(1:2,1),vertices(1:2,2),quad(:,1),quad(:,2));
        if (numel(xi) > 1)
            quadGenerated = false;
            break
        end
        i1 = 2*i + 4;
        i2 = i1 + 1;
        vertices(i1:i2,:) = quad(3:4,:);
   end
   if (~quadGenerated)
       continue
   end
   
   % check if last line segment intersects with the last quad
   xi = polyxpoly(vertices(end-1:end,1),vertices(end-1:end,2),quad(:,1),quad(:,2));
   if (numel(xi) > 1)
       %disp('last segment intersects last quad')
       continue
   end
   
   % check that last line segment doesn't intersect with any other quads
   for i = 1:(nBranches - 1)
       xi = polyxpoly(vertices(end-1:end,1),...
                      vertices(end-1:end,2),...
                      vertices(connectivity(i+1,:),1),...
                      vertices(connectivity(i+1,:),2));
       if (~isempty(xi))
           quadGenerated = false;
           break
       end
   end
   if (~quadGenerated)
       continue
   end
   generated = true;
   break
end


if (~generated && showWarning)
    warning('maximum attempts exceed, generated parallel network does not satisfy specifications')
end

end