%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Marcus on 10/29/2014
%%% Copyright 2014 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this function estimates the lengths of the arrays required for the
% assembly of the sparse stiffness matrices
% INPUT:
%   nElems: number of elements
%   nNodesPerElem: number of nodes per element
%   nIGFEMelems: number of IGFEM elements
%   estNelemsPerNode: estimated number of elements sharing a node
%   nNodesPerFace: number of nodes per face in 3D/edge in 2D
%   nDirichletNodes: number of Dirichlet nodes
function [nff,nfp,npp] = estimate_nff_nfp_npp(nElems, ...
                                              nNodesPerElem, ...
                                              nIGFEMelems, ...
                                              estNElemsPerNode, ...
                                              nNodesPerFace, ...
                                              nDirichletNodes)
nff = (nElems - nIGFEMelems)*nNodesPerElem^2 + ceil(nIGFEMelems*(2*nNodesPerElem)^2);
nfp = nDirichletNodes*estNElemsPerNode*nNodesPerFace*(nNodesPerElem-nNodesPerFace);
%npp = nDirichletNodes*estNElemsPerNode*nNodesPerFace^2;
npp = ceil(1.4*nfp);
end