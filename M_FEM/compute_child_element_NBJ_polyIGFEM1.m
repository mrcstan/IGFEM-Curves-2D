%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Marcus Tan on 8/26/2013
%%% Last modified date: 10/29/2014
%%% Copyright 2013 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this function compute the Jacobian of the child element nurbs surface 
% at a particular parameteric and the contribution of a child element to 
% the shape functions and their derivatives.
% INPUT: 
%   nurbsSurf: nurbs surface of child element
%   nurbsInd: cell array the length of the enrichment nodes of the double...
%             indices of the basis functions of the nurbs surface
%   iPON: indices of the original nodes of the parent element
%   iCENP: local number of the enrichment nodes of the child element in the parent
%          element
%   iCEN: local number of enrichment nodes of the child element wrt itself
%   nNodes: total number of nodes in parent element
%   Xel: 2x3 array of the coordinates of the original nodes of parent
%               element
%   rParams: a matrix of parametric coordinates, each col corresponding to
%            a point
%   rPhy: global coordinate
%         Provide either rParams or rPhy but not both
% OUTPUT:
%   Let nPts be the number of parametric coordinates
%   N: a matrix of shape functions with nNodes rows and nPts cols, each col corresponding
%       to a point
%   B: a matrix of shape function derivatives with Nnodes rows and 2*nPts
%       cols, every two cols corresponds to a point
%   detJ: a vector of determinants of Jacobian corresponding to each point
%   rPhy: a 2 by nPts matrix of global coordinates 
function [N,B,detJ,rPhy,calcDer] ...
                = compute_child_element_NBJ_polyIGFEM(iPON,...
                                                      iCENP,...
                                                      iCEN, ...
                                                      nNodes, ...
                                                      Xel, ...
                                                      Xch, ...
                                                      isTriangle, ...
                                                      Bel,...
                                                      calcDer,...
                                                      rPhy,...
                                                      rParams)
if ((~isempty(rParams) && ~isempty(rPhy)) || ...
    (isempty(rParams) && isempty(rPhy)))
    error('either local coordinate or global coordinate should be provided')
end
N = zeros(nNodes,1);
if (isempty(rParams))
    rParams = local_coord(rPhy,Xch);
end
if (calcDer)
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
    if (isempty(rPhy))
        rPhy = Xch*Nch;
    end
    Xloc = local_coord(rPhy,Xel);
    [N(iPON),DNel] = shape_funct(Xloc,1);
    
    if (isempty(Bel))
        Bel = DNel/(Xel*DNel);
    end
    
    B(iPON,:) = Bel;
    
    if (isTriangle)
        calcDer = false;
    end
else
    if (isTriangle)
        Nch = shape_funct(rParams,1);
    else
        Nch = shape_funct(rParams,2);
    end
    N(iCENP) = Nch(iCEN);
    rPhy = Xch*Nch;
    Xloc = local_coord(rPhy,Xel);
    N(iPON) = shape_funct(Xloc,1);
    
end




end
