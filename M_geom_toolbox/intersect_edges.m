%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Marcus Tan on 10/31/2014
%%% Last modified date: 10/31/2014
%%% Copyright 2013 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this function finds the pairwise intersections between a set of edges in
% edges1 and another set of edges in edges2
% INPUT:
%   edges1: a Nx4 matrix of N edges of the form [x11,x21,...,xN1;
%                                                y11,y21,...,yN1;
%                                                x12,x22,...,xN2;
%                                                y12,y22,...,yN2]';
%   edges2: a Nx4 matrix of N edges
%   tol: default = 1e-14;
% OUTPUT:
%   intPts: m intersection points of the form 
%           [x1,x2,...,xm;
%            y1,y2,,...,ym];
%   intPars1: intersection parametric coordinates of first set of edges
%   intPars2: intersection parametric coordinates of second set of edges
%   intInds:  vector of indices of the pairs of edges that intersects. 
%             the number of elements is m, each element corresponding to a
%             point in intPts
%   isCoincident: logical vector indicating which pairs of edges are coincident
%                 the number of elements is m, each element corresponding
%                 to a point in intPts
function [intPts, intPars1, intPars2, intInds,isCoincident]...
                = intersect_edges(edges1, edges2, tol)



% ensure input arrays are same size
N1  = size(edges1, 1);
N2  = size(edges2, 1);

% ensure input have same size
if N1~=N2
    if N1==1
        edges1 = repmat(edges1, [N2 1]);
    elseif N2==1
        edges2 = repmat(edges2, [N1 1]);
    end
end

% tolerance for precision
if (nargin < 3)
    tol = 1e-13;
end
tol1 = 1+tol;
%% Detect parallel and colinear cases

% indices of parallel edges
%par = abs(dx1.*dy2 - dx1.*dy2)<tol;
par = isParallel(edges1(:,3:4)-edges1(:,1:2), edges2(:,3:4)-edges2(:,1:2));

% indices of colinear edges
% equivalent to: |(x2-x1)*dy1 - (y2-y1)*dx1| < eps
col = find(abs(  (edges2(:,1)-edges1(:,1)).*(edges1(:,4)-edges1(:,2)) - ...
            (edges2(:,2)-edges1(:,2)).*(edges1(:,3)-edges1(:,1)) )<tol & par);

% Parallel edges have no intersection -> return [NaN NaN]
% x0(par & ~col) = NaN;
% y0(par & ~col) = NaN;


%% Process colinear edges

% colinear edges may have 0, 1 or infinite intersection
% Discrimnation based on position of edges2 vertices on edges1
nCollinears = numel(col); 
if nCollinears>0
    nCollinears2 = 2*nCollinears;
    colIntPts = nan(2,nCollinears2);
    colIntPars1 = nan(nCollinears2,1);
    colIntPars2 = nan(nCollinears2,1); 
    %intFlag = false(nCollinears,1);
    % compute position of edges1 vertices wrt edges2
    ev11e2 = edgePosition(edges1(col, 1:2), edges2(col, :));
    ev12e2 = edgePosition(edges1(col, 3:4), edges2(col, :));
    
    % compute position of edges2 vertices wrt edges1
    ev21e1 = edgePosition(edges2(col, 1:2), edges1(col, :));
    ev22e1 = edgePosition(edges2(col, 3:4), edges1(col, :));
    for i = 1:nCollinears
        intFlag = false;
        if(ev11e2(i) > -tol && ev11e2(i) < tol1);
            j = 2*i-1;
            colIntPts(:,j) = edges1(col(i),1:2)';
            colIntPars2(j) = ev11e2(i);
            colIntPars1(j) = 0;
            intFlag = true;
        end
        if(ev12e2(i) > -tol && ev12e2(i) < tol1);
            j = 2*i;
            colIntPts(:,j) = edges1(col(i),3:4)';
            colIntPars2(j) = ev12e2(i);
            colIntPars1(j) = 1;
            intFlag = true;
        end
        
        if(ev21e1(i) > -tol && ev21e1(i) < tol1);
            j = 2*i-1;
            colIntPts(:,j) = edges2(col(i),1:2)';
            colIntPars1(j) = ev21e1(i);
            colIntPars2(j) = 0;
            intFlag = true;
        end
        if(ev22e1(i) > -tol && ev22e1(i) < tol1);
            j = 2*i;
            colIntPts(:,j) = edges2(col(i),3:4)';
            colIntPars1(j) = ev22e1(i);
            colIntPars2(j) = 1;
            intFlag(i) = true;
        end
        
        if (intFlag)
            j = 2*i;
            if (colIntPars1(j-1) == colIntPars1(j) ...
                && colIntPars2(j-1) == colIntPars2(j))
                colIntPars1(j) = nan;
            end
        end
       
    end   
    %col = col(intFlag);
    del = isnan(colIntPars1);
    colIntPts(:,del) = [];
    colIntPars1(del) = [];
    colIntPars2(del) = [];
    %col = [col';col']; % changes made here
    col = [col,col]';
    col = col(:);
    col(del) = [];
    
else
    colIntPts = [];
    colIntPars1 = [];
    colIntPars2 = [];
    col = [];
end


%% Process non parallel cases

% process intersecting edges whose interecting lines intersect
nonPar = find(~par);

% use a test to avoid process empty arrays
nNonPar = numel(nonPar);
if nNonPar > 0
    % extract base parameters of supporting lines for non-parallel edges
    x1  = edges1(nonPar,1);
    y1  = edges1(nonPar,2);
    dx1 = edges1(nonPar,3)-x1;
    dy1 = edges1(nonPar,4)-y1;
    x2  = edges2(nonPar,1);
    y2  = edges2(nonPar,2);
    dx2 = edges2(nonPar,3)-x2;
    dy2 = edges2(nonPar,4)-y2;

    % compute intersection points of supporting lines
    delta = (dx2.*dy1-dx1.*dy2);
    x0 = ((y2-y1).*dx1.*dx2 + x1.*dy1.*dx2 - x2.*dy2.*dx1) ./ delta;
    y0 = ((x2-x1).*dy1.*dy2 + y1.*dx1.*dy2 - y2.*dx2.*dy1) ./ -delta;
        
    % compute position of intersection points on each edge
    % nonParIntPars1 is position on edges1, nonParIntPars2 is position on edges2
    nonParIntPars1  = ((y0-y1).*dy1 + (x0-x1).*dx1) ./ (dx1.*dx1+dy1.*dy1);
    nonParIntPars2 = ((y0-y2).*dy2 + (x0-x2).*dx2) ./ (dx2.*dx2+dy2.*dy2);

    % check position of points on edges.
    % it should be comprised between 0 and 1 for both nonParIntPars1 and nonParIntPars2.
    % if not, the edges do not intersect
    in = nonParIntPars1>-tol & nonParIntPars1<tol1 & nonParIntPars2>-tol & nonParIntPars2<tol1;
    nonParIntPts = [x0(in)';y0(in)'];
    nonPar = nonPar(in);
    nonParIntPars1 = nonParIntPars1(in);
    nonParIntPars2 = nonParIntPars2(in);
else
    nonParIntPts = [];
    nonParIntPars1 = [];
    nonParIntPars2 = [];
    nonPar = [];
end

% combine results

intPts = [colIntPts,nonParIntPts];
intPars1 = [colIntPars1;nonParIntPars1];
intPars2 = [colIntPars2;nonParIntPars2];
intInds = [col;nonPar];
isCoincident = [true(numel(col),1);false(numel(nonPar),1)];

end
