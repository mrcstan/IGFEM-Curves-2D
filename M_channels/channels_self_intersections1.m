%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Marcus Tan on 1/16/2015
%%% Copyright 2015 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this function determines if there is any self-intersection in a channel
% network, i.e, if the channels intersect at point other than the branching
% points
% REMARKS: I neglect the case where two channels of same length are 
%          exactly on top of each other with both ends coinciding
function intersectFlag = channels_self_intersections(channels,intersectTol)
if (~isfield(channels,'lineSegs'))
    warning('only implemented for linear NURBS')
    intersectFlag = false;
    return;
end
intersectTol1 = 1.0 - intersectTol;
for i = 1:(numel(channels.nurbs)-1)
    for j = 1:size(channels.lineSegs{i},1)
        [~,intPar1,intPar2] ...
            = intersect_edges(channels.lineSegs{i}(j,:),...
                              cat(1,channels.lineSegs{i+1:end}),...
                              intersectTol);
        
        if (any(intPar1 >= intersectTol & intPar1 <= intersectTol1) ...
            || any(intPar2 >= intersectTol & intPar2 <= intersectTol1))
            intersectFlag = true;
            return;
        end
        
    end
end
intersectFlag = false;
end