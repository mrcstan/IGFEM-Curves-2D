%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Marcus Tan on 4/8/2013
%%% Last modified date: 4/8/2013
%%% Copyright 2013 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate perpendicular vectors perpendicular distances from edge nodes
% to line segments
function [sNode2seg,dist_sNode2seg,eNode2seg,dist_eNode2seg]=...
      node2seg(edge_node,nodeCoord,seg,itrface)
iPt1=itrface.contvty(seg,1);
sNode=edge_node(1);     
sNode2sPt=[nodeCoord(sNode,1)-itrface.pts(iPt1,1),...
           nodeCoord(sNode,2)-itrface.pts(iPt1,2)];
proj=dot(sNode2sPt,itrface.tangent(seg,:));
projVec=proj*itrface.tangent(seg,:);
sNode2seg=projVec-sNode2sPt;
dist_sNode2seg=norm(sNode2seg,2);
eNode=edge_node(2);
eNode2sPt=[nodeCoord(eNode,1)-itrface.pts(iPt1,1),...
           nodeCoord(eNode,2)-itrface.pts(iPt1,2)];
proj=dot(eNode2sPt,itrface.tangent(seg,:));
projVec=proj*itrface.tangent(seg,:);
eNode2seg=projVec-eNode2sPt;
dist_eNode2seg=norm(eNode2seg,2);
end