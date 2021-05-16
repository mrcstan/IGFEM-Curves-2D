%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Marcus Tan on 7/24/2013
%%% Last modified date: 7/12/2014
%%% Copyright 2013 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this function determines whether a point lies in the interior of a triangle
% it is interior if it lies more than a perpendicular distance dd from the 
% edge
% Also determines if a point lies on an edge and what edge it lies on
% OUTPUT:
%   in: same length as xpt indicating if the points are in the
%               triangle interior or on the edges
%   locEdge: same length as xpt indicating which local edges the points lie on
%            If 1,2,3 are the vertices of the triangle, the local edge 1 = 
%            12, local edge 2 = 23, local edge 3 = 13. If not on an edge, locEdge = 0 
function [in,locEdge,Xloc]=inTriangle(xpt,ypt,xv,yv,dd)
npt=length(xpt);
npty=length(ypt);
if(npt~=npty)
    error('inTriangle: x-coordinates and y-coordinates must be of same length');
end
if(length(xv)~=3 || length(yv)~=3)
    error('inTriangle: must specify 3 and only 3 vertices');
end

xmin=min(xv)-dd;
xmax=max(xv)+dd;
ymin=min(yv)-dd;
ymax=max(yv)+dd;
in(1:npt,1) = false;
locEdge(1:npt,1) = 0;
Xloc = nan(2,npt);
isCandidate=(xpt>=xmin & xpt<=xmax & ypt>=ymin & ypt<=ymax);
if(any(isCandidate))
    % calculate tolerance
    triangleArea=polyarea(xv,yv);
    edge1len=norm([xv(2)-xv(1),yv(2)-yv(1)],2);
    edge2len=norm([xv(3)-xv(2),yv(3)-yv(2)],2);
    edge3len=norm([xv(1)-xv(3),yv(1)-yv(3)],2);
    factor=0.5*dd/triangleArea;
    edge1tol=factor*edge1len;
    edge2tol=factor*edge2len;
    edge3tol=factor*edge3len;
end
for i=1:npt
    if (isCandidate(i))
        %{
        v0=[xv(3)-xv(1);yv(3)-yv(1)];
        v1=[xv(2)-xv(1);yv(2)-yv(1)];
        v2=[xpt(i)-xv(1);ypt(i)-yv(1)];
        A12=v0'*v1;
        A=[v0'*v0,A12;A12,v1'*v1];
        b=[v2'*v0;v2'*v1];
        u=A\b;
        %}
        % shape=1;
        Xloc(:,i) = local_coord([xpt(i);ypt(i)],[xv';yv'],1);
        if(Xloc(1,i)>=edge3tol && Xloc(2,i)>=edge1tol ...
              && (Xloc(1,i)+Xloc(2,i))<=(1-edge2tol))
            in(i) = true;
        elseif (abs(Xloc(1,i)) < edge3tol)
            in(i) = true;
            locEdge(i) = 3;
        elseif (abs(Xloc(2,i)) < edge1tol)
            in(i) = true;
            locEdge(i) = 1;
        elseif (abs(Xloc(1,i)+Xloc(2,i)-1) < edge2tol)
            in(i) = true;
            locEdge(i) = 2;
        end
            
    end
end
end