%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Marcus Tan on 8/26/2013
%%% Last modified date: 7/15/2013
%%% Copyright 2013 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this function creates nurbs enrichment functions
function [nurbsShape,dnurbsShape]=create_nurbs_enrichment_function(nurbsSurf,nurbsInd)
nShape=length(nurbsInd);
nurbsShape(1:nShape)=nurbsSurf;
dnurbsShape(1:nShape,1:2)={nurbsSurf};
for i=1:nShape
    nBasis=length(nurbsInd{i})/2;
    for j=1:nBasis
        k=2*j-1;
        nurbsShape(i).coefs(3,nurbsInd{i}(k),nurbsInd{i}(k+1))=1;
    end
    dnurbsShape(i,:)=nrbderiv(nurbsShape(i));
end
end