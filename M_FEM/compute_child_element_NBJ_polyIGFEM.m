%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Marcus Tan on 8/26/2013
%%% Last modified date: 10/29/2014
%%% Copyright 2013 University of Illinois at Urbana-Champaign. 
%%% All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this function compute the shape functions, the derivatives of the 
% the shape functions wrt the global coordinates and the  
% Jacobian of the child element in a polynomial IGFEM element
% INPUT: 
%   iPON: indices of the original nodes of the parent element
%   iCENP: local number of the enrichment nodes of the child element in 
%          the parent element
%   iCEN: local number of enrichment nodes of the child element wrt itself
%   nNodes: total number of nodes in parent element
%   Xel: 2x3 array of the coordinates of the original nodes of parent
%               element
%   Xch: 2x3 array of the coordinates of the child element nodes
%   isTriangle: true if child element is triangular
%   Bel: matrix of the original shape function derivative, 
%   the ith row corresponds to the ith shape function and 
%   the jth col corresponds to the jth dimension
%   Xloc: local coordinates wrt child element
% OUTPUT:
%   N: shape functions
%   B: matrix of shape function derivatives, the ith row corresponds to the 
%      ith shape function and the jth col corresponds to the jth dimension
%   detJ: determinant of Jacobian
function [N,B,detJ] ...
                = compute_child_element_NBJ_polyIGFEM(iPON,...
                                                      iCENP,...
                                                      iCEN, ...
                                                      nNodes, ...
                                                      Xel, ...
                                                      Xch, ...
                                                      isTriangle, ...
                                                      Bel,...
                                                      rParams)

N = zeros(nNodes,1);
if nargout > 1
    B = zeros(nNodes,2);
    if (isTriangle)
        [Nch,DNch] = shape_funct(rParams,1);
    else
        [Nch,DNch] = shape_funct(rParams,2);
    end
    N(iCENP) = Nch(iCEN);
    Jch = Xch*DNch;
    detJ = det(Jch);
    Bch = DNch/Jch;
    B(iCENP,:) = Bch(iCEN,:);
    rPhy = Xch*Nch;
    Xloc = local_coord(rPhy,Xel,1); % parent element always assumed to be triangular
    [N(iPON),DNel] = shape_funct(Xloc,1);
    
    if (isempty(Bel))
        Bel = DNel/(Xel*DNel);
    end
    
    B(iPON,:) = Bel;
    
else
    if (isTriangle)
        Nch = shape_funct(rParams,1);
    else
        Nch = shape_funct(rParams,2);
    end
    N(iCENP) = Nch(iCEN);
    rPhy = Xch*Nch;
    Xloc = local_coord(rPhy,Xel,1); % parent element always assumed to be triangular
    N(iPON) = shape_funct(Xloc,1);    
end

end
