%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Marcus Tan on 7/13/2013
%%% Last modified date: 7/13/2013
%%% Copyright 2013 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This function finds the unit tangent of each line source segment
function tangent=find_line_source_tangent(itrface)
nTangent=size(itrface.contvty,1);
tangent=zeros(nTangent,2);
for i=1:nTangent
     tangent(i,:)=[(itrface.pts(itrface.contvty(i,2),1)...
     -itrface.pts(itrface.contvty(i,1),1)),... 
     (itrface.pts(itrface.contvty(i,2),2)...
     -itrface.pts(itrface.contvty(i,1),2))]; 
     tangent(i,:)=tangent(i,:)/norm(tangent(i,:),2);
end