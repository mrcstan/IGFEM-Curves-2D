%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Marcus Tan on 7/24/2013
%%% Last modified date: 7/24/2013
%%% Copyright 2013 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this function finds the junctions of a set of line source segments
% INPUT: 
% contvty: a nSegx2 array of the terminal points of the line source segments.
%          nSeg is the number of line source segments 
% OUTPUT: 
% junc: a length nJunc array of the terminal points which are junctions
% nSeg: number of line source segments
function [junc,nJunc,nSeg]=find_junction(contvty)
    nSeg=size(contvty,1);
    uniquePt=unique(contvty);
    nUniquePt=length(uniquePt);
    nJunc=0;
    junc=[];
    for i=1:nUniquePt
        if(nnz(contvty==uniquePt(i))>1)
            nJunc=nJunc+1;
            junc(nJunc)=uniquePt(i);
        end
    end
end