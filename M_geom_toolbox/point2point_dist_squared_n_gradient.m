%%% Created by Marcus Tan on 1/18/2015
%%% Copyright 2015 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this function calculates the square of the pairwise distance between 
% two sets of points and the gradient with respect to each of the
% coordinate of the points
% INPUT:
%   pts1 = [x11,y11; x12,y12; ...;x1n,y1n];
%   pts2 = [x21,y21; x22,y22; ...;x2n,y2n];
% OUTPUT:
%   distSq: vector of length n, its ith entry corresponding the square of
%           the distance between point(x1i,y1i) and point (x2i,y2i)
%   distSqGrad: n x 4 matrix, i1 entry corresponding to derivative of
%               distance square wrt x1i, i2 entry "" wrt y1i, 
%               i3 entry "" wrt x2i and i4 entry "" wrt y2i 
function [distSq,distSqGrad] = point2point_dist_squared_n_gradient(pts1,pts2)
if (size(pts1,1) ~= size(pts2,1))
    error('pts1 and pts2 must have same number of rows')
end
if (size(pts1,2) ~= 2 || size(pts2,2) ~= 2)
    error('pts1 and pts2 must have two columns')
end
distSq = (pts1(:,1) - pts2(:,1)).^2 + (pts1(:,2) - pts2(:,2)).^2;
if (nargout > 1)
    distSqGrad = [2*(pts1(:,1) - pts2(:,1)),2*(pts1(:,2) - pts2(:,2))...
                  2*(pts2(:,1) - pts1(:,1)),2*(pts2(:,2) - pts1(:,2))];
end
end