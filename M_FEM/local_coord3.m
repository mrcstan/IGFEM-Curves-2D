%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Marcus Tan on 8/21/2014
%%% Copyright 2012 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%This function returns barycentric coordinates (except the last dependent
%one) for a tetrahedral or a triangular element (in 2D or 3D)
% INPUT:
% Xglo: a matrix of the global coordinates. each col corresponds to a point
% Xel a 3x4 matrix of the coordinates of the vertices of the tetrahedral
%      element
%      or 2x3 or 3x3 matrix of the coordinates of the vertices a triangle
%      in 2D or 3D respectively.
%      Note that Xglo must have the same number of rows as Xel
function Xloc = local_coord(Xglo, Xel)
if (size(Xel,1) ~= size(Xglo,1))
    error('element dimension must be the same as global coordinate dimension')
end
if (size(Xel,2) == 4)
    A = [Xel(:,2)-Xel(:,1),Xel(:,3)-Xel(:,1),Xel(:,4)-Xel(:,1)];
elseif (size(Xel,2) == 3)
    A = [Xel(:,2)-Xel(:,1),Xel(:,3)-Xel(:,1)];
else
    error('Input 2 must either have 3 or 4 columns')
end
rhs = bsxfun(@minus,Xglo,Xel(:,1));
%fprintf('condition number = %g \n', cond(A));
Xloc=A\rhs; 
end

