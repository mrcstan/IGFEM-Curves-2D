%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Marcus Tan on 12/2/2014
%%% Last modified date: 12/2/2014
%%% Copyright 2014 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this function calculates the area and the derivative wrt to the locations
% of the polygon vertices
function [area,dArea] = polygon_area(X,Y)
%{
area = polyarea(X,Y);
if (ispolycw(X,Y))
    area = -area;
end
%}
area = polygonArea([X(:),Y(:)]);
if (nargout > 1)
    nv = numel(X);
    dArea = nan(2*nv,1);
    ind = [nv,1:nv,1]; % wrapped index
    for i = 2:nv+1
        dArea(2*ind(i)-1) = 0.5*(Y(ind(i+1))-Y(ind(i-1)));
        dArea(2*ind(i)) = 0.5*(X(ind(i-1))-X(ind(i+1)));
    end
   
end
end