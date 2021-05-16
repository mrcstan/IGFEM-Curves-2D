%%% Created by Marcus Tan on 1/18/2015
%%% Copyright 2015 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this function calculates the square of the pairwise distance between 
% two sets of points and the gradient with respect to each of the
% coordinate of the points
% INPUT:
%   points = [x1,y1; x2,y2; ...;xn,yn];
%   edges = [x11,y11,x12,y12; x21,y21,x22,y22; ...;xn1,yn1,xn2,yn2];
% OUTPUT:
%   distSq: vector of length n, its ith entry corresponding the square of
%           the distance between point(x1i,y1i) and point (x2i,y2i)
%   distSqGrad: n x 4 matrix, i1 entry corresponding to derivative of
%               distance square wrt x1i, i2 entry "" wrt y1i, 
%               i3 entry "" wrt x2i and i4 entry "" wrt y2i 
function [distSq,distSqGrad] = point2edge_dist_squared_n_gradient(points, edges)
% direction vector of each edge
dx = edges(:, 3) - edges(:,1);
dy = edges(:, 4) - edges(:,2);

% compute position of points projected on the supporting line
% (Size of tp is the max number of edges or points)
delta = dx .* dx + dy .* dy;
tp = ((points(:, 1) - edges(:, 1)) .* dx + (points(:, 2) - edges(:, 2)) .* dy) ./ delta;

% ensure degenerated edges are correclty processed (consider the first
% vertex is the closest)
tp(delta < eps) = 0;

% change position to ensure projected point is located on the edge
tp(tp < 0) = 0;
tp(tp > 1) = 1;

% coordinates of projected point
p0 = [edges(:,1) + tp .* dx, edges(:,2) + tp .* dy];

% compute distance between point and its projection on the edge
distSq = (points(:,1) - p0(:,1)) .^ 2 + (points(:,2) - p0(:,2)) .^ 2;

if nargout > 1
    
end
