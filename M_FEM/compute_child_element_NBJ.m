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
%   iCENP: local number of the enrichment node of the child element in the parent
%          element
%   nNodes: total number of nodes in parent element
%   Xel: 2x3 array of the coordinates of the original nodes of parent
%               element
%   rParams: a matrix of parametric coordinates, each col corresponding to
%            a point
% OUTPUT:
%   Let nPts be the number of parametric coordinates
%   N: a matrix of shape functions with nNodes rows and nPts cols, each col corresponding
%       to a point
%   B: a matrix of shape function derivatives with Nnodes rows and 2*nPts
%       cols, every two cols corresponds to a point
%   detJ: a vector of determinants of Jacobian corresponding to each point
%   rPhy: a 2 by nPts matrix of global coordinates 
function [N,B,detJ,rPhy]=compute_child_element_NBJ(nurbsSurf,...
                                                   dnurbsSurf,...
                                                   nurbsInd,...
                                                   iPON,...
                                                   iCENP,...
                                                   nNodes, ...
                                                   Xel,...
                                                   Bel,...
                                                   calcDer,...
                                                   rParams)
                                               
nPts = size(rParams,2);
N = zeros(nNodes,nPts);


% nurbs basis functions
basis = nrbsrfbasisfun(rParams,nurbsSurf);
nInds = numel(nurbsInd);

% calculate B
if (calcDer)
    % global coordinates and Jacobian
    %rPhy = nrbeval(nurbsSurf,rParams);
    [rPhy,Jp2p] = nrbdeval(nurbsSurf,dnurbsSurf,rParams); 
    rPhy = rPhy(1:2,:);

    % local coordinates in original element
    rLoc = local_coord(rPhy, Xel,1); % assume original element is triangular

    % contribution to original nodes
    if (isempty(Bel))
        [Nel, DNel] = shape_funct(rLoc, 1);
        Bel = DNel/(Xel*DNel);
    else
        Nel = shape_funct(rLoc, 1);
    end
    N(iPON,:) = Nel;
    
    detJ = zeros(nPts,1); 
    DNch = zeros(nInds,2*nPts);
    [dBasis_xi,dBasis_eta] = nrbsrfbasisfunder(rParams,nurbsSurf);
    for i = 1:nInds
        N(iCENP(i),:) = sum(basis(:,nurbsInd{i}),2)';
        DNch(i,:) = reshape([sum(dBasis_xi(:,nurbsInd{i}),2)';...
                            sum(dBasis_eta(:,nurbsInd{i}),2)'],1,[]);
    end
    B = zeros(nNodes,2*nPts);
    B(iPON,:) = repmat(Bel,[1,nPts]);
    J = reshape([Jp2p{1}(1:2,:);Jp2p{2}(1:2,:)],2,[]);

    for i = 1:nPts
        j = 2*i - 1;
        detJ(i) = det(J(:,j:j+1));
        %Jinv = 1/det(J)*[J(2,2),-J(1,2);-J(2,1),J(1,1)];
        %B(iCENP,j:j+1) = DNch(:,j:j+1)*Jinv;
        B(iCENP,j:j+1) = DNch(:,j:j+1)/J(:,j:j+1);
    end

else
    % global coordinates
    rPhy = nrbeval(nurbsSurf,rParams);
    rPhy = rPhy(1:2,:);

    % local coordinates in original element
    rLoc = local_coord(rPhy, Xel, 1);

    N(iPON,:) = shape_funct(rLoc, 1);
    for i = 1:nInds
        N(iCENP(i),:) = sum(basis(:,nurbsInd{i}),2)';
    end
    B = [];
end
end
