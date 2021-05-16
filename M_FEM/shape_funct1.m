%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Ahmad Raeisi Najafi on 10/08/2012
%%% Modified by Marcus Tan
%%% Last modified date: 07/15/2013
%%% Copyright 2012 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT:
% shape: 1 - triangular element
%        2 - quadrilateral element
function [N, DN] = shape_funct(r, shape)

%This function returns the values of the triangular 2D shape functions and
% their gradient, for a point r in the reference element.

if shape == 1
    
    % This function returns the values of the triangle 2D shape functions and
    % their gradient, for a point r in the reference element.
    
    % form shape function matrix
    N = [1 - r(1) - r(2);...
         r(1); ...
         r(2)];
    
    % form shape function gradient matrix
    DN = [-1, -1; ...
          1, 0; ...
          0, 1];
elseif shape == 2

    % This function returns the values of the bilinear 2D shape functions and
    % their gradient, for a point r in the reference element.

    % form shape function matrix
    N = [1/4*(1 - r(1))*(1 - r(2)), ....
        1/4*(1 + r(1))*(1 - r(2)), ....
        1/4*(1 + r(1))*(1 + r(2)), ....
        1/4*(1 - r(1))*(1 + r(2))]';

    % form shape function gradient matrix
    DN = [ -1/4*(1 - r(2)), -1/4*(1 - r(1)); ...
        1/4*(1 - r(2)), -1/4*(1 + r(1)); ...
        1/4*(1 + r(2)),  1/4*(1 + r(1)); ...
        -1/4*(1 + r(2)),  1/4*(1 - r(1))];
end

return