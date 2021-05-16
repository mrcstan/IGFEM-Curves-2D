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
% shape: 1: triangle
%        2: quadrilateral
function Xloc = local_coord(Xglo, Xel, shape)
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
elseif (shape == 2)
    %
    Xloc = nan(size(Xglo));
    options  =  optimset( 'Algorithm', 'trust-region-dogleg',...  % Choose the algorithm to solve the nonlinear function f(x) - ax - b = 0
                         'TolX',2e-16,...      % tolerance for intersect point X_intrsect 
                         'TolFun', 2e-16,'Display','off');     % tolerence on the nonlinear funciton we solve
    for i = 1:size(Xglo,2)
        func = @(r) Xglo - Xel*shape_funct(r, shape);                 
        [Xloc(:,i), ~, exitflag] = fsolve(func,[0.5;0.5],options);
        if(exitflag<1)
            warning('exit flag < 1 for fsolve')
        end
    end
    %
    %{
    tol = 1e-10;
    tol1 = -1-tol;
    tol2 = 1+tol;
    Xloc = nan(size(Xglo));
    for i = 1:size(Xglo,2)
        % solution given by "An inverse transformation for quadrilateral
        % isoparametric elements" by Chongyu Hua
        a = Xel*[1,-1,1,-1]';
        b = Xel*[-1,1,1,-1]';
        c = Xel*[-1,-1,1,1]';
        d = 4*Xglo(:,i) - sum(Xel,2);
        ab = det([a,b]);
        ac = det([a,c]);
        ad = det([a,d]);
        bc = det([b,c]);
        bd = det([b,d]);
        cd = det([c,d]);
        
        if (   abs(a(1)*a(2)*ab*ac) > tol ...
            || (abs(a(1)) < tol && abs(a(2)*c(1)) > tol) ...
            || (abs(a(2)) < tol && abs(a(1)*b(2)) > tol))
            rts = roots([ab,-(bc+ad),-cd]);
            if (rts(1) >= tol1 && rts(1) <= tol2)
                Xloc(1,i) = rts(1);
            elseif (rts(2) >= tol1 && rts(2) <= tol2)
                Xloc(1,i) = rts(2);
            end
            Xloc(2,i) = (ad - ab*Xloc(1,i))/ac;
        elseif (abs(a(1)*a(2)) > tol && abs(ab) < tol)
            Xloc(1,i) = -a(1)*cd/(b(1)*ac+a(1)*ad);
            Xloc(2,i) = ad/ac;
        elseif (abs(a(1)*a(2)) > tol && abs(ac) < tol)
            Xloc(1,i) = ad/ab;
            Xloc(2,i) = -a(1)*bd/(c(1)*ab + a(1)*ad);
        else
            Xloc(1,i) = -cd/(a(1)*d(2) + bc);
            Xloc(2,i) = bd/(a(2)*d(1) + bc);
        end
    end
    %}
elseif (shape == -2) % quadrilateral and global coordinate is along one of the edges
    Xloc = nan(size(Xglo));
    tol = 1e-13;
    for i = 1:size(Xglo,2)
        % guess: it is along edge 1, eta = -1
        xi1 = (2*Xglo(1) - Xel(1,1) - Xel(1,2))/(Xel(1,2) - Xel(1,1));
        xi2 = (2*Xglo(2) - Xel(2,1) - Xel(2,2))/(Xel(2,2) - Xel(2,1));
        if (abs(xi1-xi2) < tol)
            Xloc(1,i) = 0.5*(xi1+xi2);
            Xloc(2,i) = -1;
        end
        % guess: it is along edge 2, xi = 1
        eta1 = (2*Xglo(1) - Xel(1,2) - Xel(1,3))/(Xel(1,3) - Xel(1,2));
        eta2 = (2*Xglo(2) - Xel(2,2) - Xel(2,3))/(Xel(2,3) - Xel(2,2));
        if (abs(eta1-eta2) < tol)
            Xloc(1,i) = 1;
            Xloc(2,i) = 0.5*(eta1+eta2);
        end
        % guess: it is along edge 3, eta = 1
        xi1 = (2*Xglo(1) - Xel(1,3) - Xel(1,4))/(Xel(1,3) - Xel(1,4));
        xi2 = (2*Xglo(2) - Xel(2,3) - Xel(2,4))/(Xel(2,3) - Xel(2,4));
        if (abs(xi1-xi2) < tol)
            Xloc(1,i) = 0.5*(xi1+xi2);
            Xloc(2,i) = 1;
        end
        % guess: it is along edge 4, xi = -1
        eta1 = (2*Xglo(1) - Xel(1,1) - Xel(1,4))/(Xel(1,4) - Xel(1,1));
        eta2 = (2*Xglo(2) - Xel(2,1) - Xel(2,4))/(Xel(2,4) - Xel(2,1));
        if (abs(eta1-eta2) < tol)
            Xloc(1,i) = -1;
            Xloc(2,i) = 0.5*(eta1+eta2);
        end
        warning('global coordinate does not lie on the edge of the quadrilateral element')
    end
else
    error('unknown shape')
end
end

