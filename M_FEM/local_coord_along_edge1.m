%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Marcus Tan on 12/27/2014
%%% Copyright 2012 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%This function returns barycentric coordinates (except the last dependent
%one) for a tetrahedral or a triangular element (in 2D or 3D) or
%quadrilateral element in 2D when the global coordinate is known to be
%along an edge
% INPUT:
% Xglo: a matrix of the global coordinates. each col corresponds to a point
% Xel a 3x4 matrix of the coordinates of the vertices of the tetrahedral
%      element
%      or 2x3 or 3x3 matrix of the coordinates of the vertices a triangle
%      in 2D or 3D respectively.
%      Note that Xglo must have the same number of rows as Xel
% shape: 1: triangle
%        2: quadrilateral
function Xloc = local_coord_along_edge(Xglo, Xel, shape)
dim = size(Xglo,1);
if (dim ~= 2 && dim ~= 3)
    error('must have either dim 2 or 3, i.e., 2 or 3 rows')
end
if (size(Xel,1) ~= dim)
    error('element dimension must be the same as global coordinate dimension')
end
if (shape == 1)
    if (size(Xel,2) == 4 && dim == 3)
        A = [Xel(:,2)-Xel(:,1),Xel(:,3)-Xel(:,1),Xel(:,4)-Xel(:,1)];
    elseif (size(Xel,2) == 3 && (dim == 2 || dim == 1))
        A = [Xel(:,2)-Xel(:,1),Xel(:,3)-Xel(:,1)];
    else
        error('Xel does not have the proper number of cols according to the dimension of Xglo')
    end
    rhs = bsxfun(@minus,Xglo,Xel(:,1));
    %fprintf('condition number = %g \n', cond(A));
    Xloc=A\rhs;
elseif (shape == 2) % quadrilateral and global coordinate is along one of the edges
    Xloc = nan(size(Xglo));
    tol = 1e-10;
    for i = 1:size(Xglo,2)
        % guess: it is along edge 1, eta = -1
        xi1 = (2*Xglo(1,i) - Xel(1,1) - Xel(1,2))/(Xel(1,2) - Xel(1,1));
        xi2 = (2*Xglo(2,i) - Xel(2,1) - Xel(2,2))/(Xel(2,2) - Xel(2,1));
        if (abs(xi1-xi2) < tol)
            Xloc(1,i) = 0.5*(xi1+xi2);
            Xloc(2,i) = -1;
            continue
        end
        % guess: it is along edge 2, xi = 1
        eta1 = (2*Xglo(1,i) - Xel(1,2) - Xel(1,3))/(Xel(1,3) - Xel(1,2));
        eta2 = (2*Xglo(2,i) - Xel(2,2) - Xel(2,3))/(Xel(2,3) - Xel(2,2));
        if (abs(eta1-eta2) < tol)
            Xloc(1,i) = 1;
            Xloc(2,i) = 0.5*(eta1+eta2);
            continue
        end
        % guess: it is along edge 3, eta = 1
        xi1 = (2*Xglo(1,i) - Xel(1,3) - Xel(1,4))/(Xel(1,3) - Xel(1,4));
        xi2 = (2*Xglo(2,i) - Xel(2,3) - Xel(2,4))/(Xel(2,3) - Xel(2,4));
        if (abs(xi1-xi2) < tol)
            Xloc(1,i) = 0.5*(xi1+xi2);
            Xloc(2,i) = 1;
            continue
        end
        % guess: it is along edge 4, xi = -1
        eta1 = (2*Xglo(1,i) - Xel(1,1) - Xel(1,4))/(Xel(1,4) - Xel(1,1));
        eta2 = (2*Xglo(2,i) - Xel(2,1) - Xel(2,4))/(Xel(2,4) - Xel(2,1));
        if (abs(eta1-eta2) < tol)
            Xloc(1,i) = -1;
            Xloc(2,i) = 0.5*(eta1+eta2);
            continue
        end
        warning('global coordinate does not lie on the edge of the quadrilateral element')
    end
else
    error('unknown shape')
end
end

