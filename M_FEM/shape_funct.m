%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Ahmad Raeisi Najafi on 10/08/2012
%%% Modified by Marcus Tan
%%% Last modified date: 10/28/2013
%%% Copyright 2012 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This function returns the values of the triangular 2D shape functions and
% their gradient, for a point r in the reference element
% INPUT:
% shape: 1 - triangular element
%        2 - quadrilateral element
% r: a matrix of local coordinates, each column corresponds to a point
% OUTPUT:
%   N: a matrix of shape functions, each column corresponding to a point
%   DN: if shape == 1, the matrix of shape function derivatives evaluated
%       at just one point is output
%       if shape == 2, consecutive blocks of 4x2 of the matrix correspond 
%                      to each point
function [N, DN] = shape_funct(r, shape)

if (size(r,1) ~= 2) 
    error('invalid form for local coordinates')
end
if shape == 1
    % This function returns the values of the triangle 2D shape functions and
    % their gradient, for a point r in the reference element.      
    N = [1 - r(1,:) - r(2,:); ...
         r(1,:);
         r(2,:)];
    % form shape function gradient matrix
    if (nargout > 1)
        DN = [-1, -1; ...
              1, 0; ...
              0, 1];
        % DN = repmat(DN,[1,size(r,2)]);
    end
elseif shape == 2

    % This function returns the values of the bilinear 2D shape functions and
    % their gradient, for a point r in the reference element.

    % form shape function matrix
    %{
    N = [1/4*(1 - r(1))*(1 - r(2)), ....
        1/4*(1 + r(1))*(1 - r(2)), ....
        1/4*(1 + r(1))*(1 + r(2)), ....
        1/4*(1 - r(1))*(1 + r(2))]';
    %}
    N = 0.25 * [(1 - r(1,:)).*(1 - r(2,:));
                (1 + r(1,:)).*(1 - r(2,:));
                (1 + r(1,:)).*(1 + r(2,:));
                (1 - r(1,:)).*(1 + r(2,:))];
    % form shape function gradient matrix
    %{
    DN = [ -1/4*(1 - r(2)), -1/4*(1 - r(1)); ...
        1/4*(1 - r(2)), -1/4*(1 + r(1)); ...
        1/4*(1 + r(2)),  1/4*(1 + r(1)); ...
        -1/4*(1 + r(2)),  1/4*(1 - r(1))];
    %}
    if (nargout > 1)
        DN = 0.25*[reshape([(r(2,:) - 1);(r(1,:) - 1)], 1, []); 
                   reshape([(1 - r(2,:));-(r(1,:) + 1)], 1, []);
                   reshape([(r(2,:) + 1);(r(1,:) + 1)], 1, []);
                   reshape([-(r(2,:) + 1);(1-r(1,:))], 1, [])];
    end
end

end