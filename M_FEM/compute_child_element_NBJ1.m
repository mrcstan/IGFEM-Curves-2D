%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Marcus Tan on 8/26/2013
%%% Last modified date: 7/17/2013
%%% Copyright 2013 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this function compute the Jacobian of the child element nurbs surface 
% at a particular parameteric and the contribution of a child element to 
% the shape functions and their derivatives.
% INPUT: 
%   nurbsSurf: nurbs surface of child element
%   dnurbsSurf: derivative of nurbs surface of child element  
%   nurbsShapeInd: a 2xn pairs indices of the basis functions of the nurbs surface
%                  which are enrichment functions of the enrichment nodes
%   iPON: indices of the original nodes of the parent element
%   iCENP: local number of the enrichment node of the child element in the parent
%          element
%   nNode: total number of nodes in parent element
%   xOriginal: 2x3 array of the coordinates of the original nodes of parent
%               element
%   rParam: parametric coordinate where the outputs are to be evaluated
% OUTPUT:
%   detJ:
%   N:
%   B:
%   rPhy:
function [N,B,detJ,rPhy]=compute_child_element_NBJ(nurbsSurf,dnurbsSurf,...
                                         nurbsShape,dnurbsShape,iPON,iCENP,...
                                         nNode,xOriginal,rParam)

N=zeros(nNode,1);
B=zeros(nNode,2);

% contribution to parent element shape functions
rPhy=nrbeval(nurbsSurf,[rParam(1);rParam(2)]);
rPhy=rPhy(1:2);
original_shape=1;
rLoc = local_coord(original_shape, rPhy, xOriginal);
[N(iPON), DN] = shape_funct(rLoc, original_shape);
Joriginal=xOriginal*DN;
B(iPON,1:2)=DN/Joriginal;

% contribution of child element to shape function
[~,Jp2p]=nrbdeval(nurbsSurf,dnurbsSurf,[rParam(1);rParam(2)]);
Jp2p=cell2mat(Jp2p);
Jp2p=Jp2p(1:2,1:2);


[N(iCENP), DN] = nurbs_shape_func_eval(nurbsShape,dnurbsShape,rParam);
if (nargout > 1 )    
    B(iCENP,1:2)=DN/Jp2p;
    detJ=det(Jp2p);
end
end
