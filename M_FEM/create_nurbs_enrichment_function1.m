%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Marcus Tan on 8/26/2013
%%% Last modified date: 9/19/2013
%%% Copyright 2013 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this function creates nurbs enrichment functions
function [nurbsShape,dnurbsShape]=create_nurbs_enrichment_function(nurbsSurf,nurbsInd,...
                                                                specialInd)
nShape=length(nurbsInd);
nurbsShape(1:nShape)=nurbsSurf;
dnurbsShape(1:nShape,1:2)={nurbsSurf};
for i=1:nShape
    nBasis=length(nurbsInd{i})/2;
    for j=1:nBasis
        k=2*j-1;
        nurbsShape(i).coefs(3,nurbsInd{i}(k),nurbsInd{i}(k+1))=1;
        if(nargin>2 && specialInd(i))
            if(nurbsInd{i}(k)==1 && nurbsSurf.number(2)>2)
                val=linspace(0,1,nurbsSurf.number(2));
                val=val(2:end-1);
                nurbsShape(i).coefs(3,1,2:end-1)=val;
            elseif(nurbsInd{i}(k+1)==nurbsSurf.number(2) && nurbsSurf.number(1)>2)
                val=linspace(0,1,nurbsSurf.number(1));
                val=val(2:end-1);
                nurbsShape(i).coefs(3,2:end-1,nurbsSurf.number(2))=val;
            end
        end
    end
    dnurbsShape(i,:)=nrbderiv(nurbsShape(i));
end
end