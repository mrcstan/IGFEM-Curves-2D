function indicatorMatrix = line_segment_intersect_indicator(XY1,XY2,tol)
%LINESEGMENTINTERSECT Intersections of line segments.
%   OUT = LINESEGMENTINTERSECT(XY1,XY2) finds the 2D Cartesian Coordinates of
%   intersection points between the set of line segments given in XY1 and XY2.
%
%   XY1 and XY2 are N1x4 and N2x4 matrices. Rows correspond to line segments. 
%   Each row is of the form [x1 y1 x2 y2] where (x1,y1) is the start point and 
%   (x2,y2) is the end point of a line segment:
%
%                  Line Segment
%       o--------------------------------o
%       ^                                ^
%    (x1,y1)                          (x2,y2)
%
%   OUT is a structure with fields:
%
%   indicatorMatrix : N1xN2 indicator matrix where the entry (i,j) is 1 if
%       line segments XY1(i,:) and XY2(j,:) intersect.

% Version: 1.00, April 03, 2010
% Version: 1.10, April 10, 2010
% Author:  U. Murat Erdem

% CHANGELOG:
%
% Ver. 1.00: 
%   -Initial release.
% 
% Ver. 1.10:
%   - Changed the input parameters. Now the function accepts two sets of line
%   segments. The intersection analysis is done between these sets and not in
%   the same set.
%   - Changed and added fields of the output. Now the analysis provides more
%   information about the intersections and line segments.
%   - Performance tweaks.

% I opted not to call this 'curve intersect' because it would be misleading
% unless you accept that curves are pairwise linear constructs.
% I tried to put emphasis on speed by vectorizing the code as much as possible.
% There should still be enough room to optimize the code but I left those out
% for the sake of clarity.
% The math behind is given in:
%   http://local.wasp.uwa.edu.au/~pbourke/geometry/lineline2d/
% If you really are interested in squeezing as much horse power as possible out
% of this code I would advise to remove the argument checks and tweak the
% creation of the OUT a little bit.

%%% Argument check.
% Modified by Marcus Tan to output only the indicator matrix
%-------------------------------------------------------------------------------

validateattributes(XY1,{'numeric'},{'2d','finite'});
validateattributes(XY2,{'numeric'},{'2d','finite'});

[n_rows_1,n_cols_1] = size(XY1);
[n_rows_2,n_cols_2] = size(XY2);

if n_cols_1 ~= 4 || n_cols_2 ~= 4
    error('Arguments must be a Nx4 matrices.');
end

if (nargin < 3)
    tol = 1e-13;
end
tol1 = 1.0 + tol;
%%% Prepare matrices for vectorized computation of line intersection points.
%-------------------------------------------------------------------------------
X1 = repmat(XY1(:,1),1,n_rows_2);
X2 = repmat(XY1(:,3),1,n_rows_2);
Y1 = repmat(XY1(:,2),1,n_rows_2);
Y2 = repmat(XY1(:,4),1,n_rows_2);

XY2 = XY2';

X3 = repmat(XY2(1,:),n_rows_1,1);
X4 = repmat(XY2(3,:),n_rows_1,1);
Y3 = repmat(XY2(2,:),n_rows_1,1);
Y4 = repmat(XY2(4,:),n_rows_1,1);

X4_X3 = (X4-X3);
Y1_Y3 = (Y1-Y3);
Y4_Y3 = (Y4-Y3);
X1_X3 = (X1-X3);
X2_X1 = (X2-X1);
Y2_Y1 = (Y2-Y1);

numerator_a = X4_X3 .* Y1_Y3 - Y4_Y3 .* X1_X3;
numerator_b = X2_X1 .* Y1_Y3 - Y2_Y1 .* X1_X3;
denominator = Y4_Y3 .* X2_X1 - X4_X3 .* Y2_Y1;

u_a = numerator_a ./ denominator;
u_b = numerator_b ./ denominator;

% Find the adjacency matrix A of intersecting lines.
%INT_X = X1+X2_X1.*u_a;
%INT_Y = Y1+Y2_Y1.*u_a;
INT_B = (u_a >= -tol) & (u_a <= tol1) & (u_b >= -tol) & (u_b <= tol1);
%PAR_B = abs(denominator) < tol;
COINC_B = (abs(numerator_a) < tol & ...
           abs(numerator_b) < tol & ...
           abs(denominator) < tol);


% Arrange output.
%{
out.intAdjacencyMatrix = INT_B;
out.intMatrixX = INT_X .* INT_B;
out.intMatrixY = INT_Y .* INT_B;
out.intNormalizedDistance1To2 = u_a;
out.intNormalizedDistance2To1 = u_b;
out.parAdjacencyMatrix = PAR_B;
out.coincAdjacencyMatrix= COINC_B;
%}

indicatorMatrix = INT_B | COINC_B;
end

