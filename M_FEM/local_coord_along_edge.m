%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Marcus Tan on 12/27/2014
%%% Copyright 2012 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%This function returns barycentric coordinates or a triangular element or
%quadrilateral element in 2D when the global coordinate is known to be
%along an edge
% INPUT:
% Xglo: a matrix of the global coordinates. each col corresponds to a point
% Xel a 2x3 matrix of the coordinates of the vertices of a triangle
%      or 2x4 matrix of the coordinates of the vertices of a quadrilateral
%      Note that Xglo must have the same number of rows as Xel
% shape: 1: triangle
%        2: quadrilateral
function Xloc = local_coord_along_edge(Xglo, Xel, shape, edge)
dim = size(Xglo,1);
if (dim ~= 2 && dim ~= 3)
    error('must have either dim 2 or 3, i.e., 2 or 3 rows')
end
if (size(Xel,1) ~= dim)
    error('element dimension must be the same as global coordinate dimension')
end
tol = 1e-13;
Xloc = nan(size(Xglo));
if (shape == 1) % triangular and global coordinate is along one of the edges
    for i = 1:size(Xglo,2)
        switch edge
            case 1 % eta = 0
                denom1 = Xel(1,2)-Xel(1,1);
                denom2 = Xel(2,2)-Xel(2,1);       
                if (abs(denom1) > abs(denom2))
                    Xloc(1,i) = (Xglo(1,i)-Xel(1,1))/denom1;
                    Xloc(2,i) = 0;
                else
                    Xloc(1,i) = (Xglo(2,i)-Xel(2,1))/denom2;
                    Xloc(2,i) = 0;
                end
            case 2 % xi+eta = 1
                denom1 = Xel(1,2)-Xel(1,3);
                denom2 = Xel(2,2)-Xel(2,3);
                if (abs(denom1) > abs(denom2))
                    Xloc(1,i) = (Xglo(1,i)-Xel(1,3))/denom1;
                    Xloc(2,i) = 1-Xloc(1,i);
                else
                    Xloc(1,i) = (Xglo(2,i)-Xel(2,3))/denom2;
                    Xloc(2,i) = 1-Xloc(1,i);
                end    
            case 3 % xi =0
                denom1 = Xel(1,3)-Xel(1,1);
                denom2 = Xel(2,3)-Xel(2,1);
                if (abs(denom1) > abs(denom2))
                    Xloc(1,i) = 0;
                    Xloc(2,i) = (Xglo(1,i)-Xel(1,1))/denom1;
                else
                    Xloc(1,i) = 0;
                    Xloc(2,i) = (Xglo(2,i)-Xel(2,1))/denom2;
                end          
            otherwise
                error('triangle: unknown edge')
                
        end
        if (abs(denom1) < tol && abs(denom2) < tol)
             warning('small denominators %g, %g for triangle edge %i',denom1,denom2,edge)
        end
    end   
elseif (shape == 2) % quadrilateral and global coordinate is along one of the edges
        
    for i = 1:size(Xglo,2)
        switch edge
            case 1
                denom1 = Xel(1,2) - Xel(1,1);
                denom2 = Xel(2,2) - Xel(2,1);
                if (abs(denom1) > abs(denom2))
                    Xloc(1,i) = (2*Xglo(1,i) - Xel(1,1) - Xel(1,2))/denom1;
                    Xloc(2,i) = -1;
                else
                    Xloc(1,i) = (2*Xglo(2,i) - Xel(2,1) - Xel(2,2))/denom2;
                    Xloc(2,i) = -1;
                end
            case 2
                denom1 = Xel(1,3) - Xel(1,2);
                denom2 = Xel(2,3) - Xel(2,2);
                if (abs(denom1) > abs(denom2))
                    Xloc(1,i) = 1;
                    Xloc(2,i) = (2*Xglo(1,i) - Xel(1,2) - Xel(1,3))/denom1;
                else
                    Xloc(1,i) = 1;
                    Xloc(2,i) = (2*Xglo(2,i) - Xel(2,2) - Xel(2,3))/denom2;
                end
            case 3
                denom1 = Xel(1,3) - Xel(1,4);
                denom2 = Xel(2,3) - Xel(2,4);
                if (abs(denom1) > abs(denom2))
                    Xloc(1,i) = (2*Xglo(1,i) - Xel(1,3) - Xel(1,4))/denom1;
                    Xloc(2,i) = 1;
                else
                    Xloc(1,i) = (2*Xglo(2,i) - Xel(2,3) - Xel(2,4))/denom2;
                    Xloc(2,i) = 1;
                end
            case 4
                denom1 = Xel(1,4) - Xel(1,1);
                denom2 = Xel(2,4) - Xel(2,1);
                if (abs(denom1) > abs(denom2))
                    Xloc(1,i) = -1;
                    Xloc(2,i) = (2*Xglo(1,i) - Xel(1,1) - Xel(1,4))/denom1;
                else
                    Xloc(1,i) = -1;
                    Xloc(2,i) =  (2*Xglo(2,i) - Xel(2,1) - Xel(2,4))/denom2;
                end
            otherwise
                error('unknown edge')
            
        end
        if (abs(denom1) < tol && abs(denom2) < tol)
             warning('small denominators %g, %g for for quad edge %i',denom1,denom2,edge)
        end
    end
else
    error('unknown shape')
end
end

