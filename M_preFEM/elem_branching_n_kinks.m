%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Marcus Tan on 3/22/2016
%%% Copyright 2013 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this function determine the branching points kinks of a 
% channel network that are located in the interior or edge of each element. 
% Kinks are only present for channels consisting of straight 
% line segments
% OUTPUT:
%       elemKinks(i).kinkNums: kinks in element interior
%       elemKinks(i).nodes: nodes asscoiated with kinks in element interior
%       nodeCoords: updated nodal coordinates
%       nNodes: new total number of nodes
%       itrsect: updated face intersection info


function [elem_branch_kinks, ...
          nodeCoords, ...
          nNodes, ...
          itrsect] = elem_branching_n_kinks(DT,...
                                            elem_edges,...
                                            itrsect, ...
                                            nodeCoords, ...
                                            nNodes, ...
                                            branch_kinks, ...
                                            chanItrsectParams,...
                                            tolParams, ...
                                            tolBary,...
                                            chanDsgnParamNums, ...
                                            calcItrsectVel)
nElems = size(DT.ConnectivityList,1);                                        
elem_branch_kinks = struct('bkNums',cell(nElems,1),'nodes',cell(nElems,1));

nBKs = size(branch_kinks.XX,2);
if nBKs < 1
    return
end

% compare channel intersection parameters with those of
% the branch point or kink parameters to determine
% if they have been found as intersections and remove them from the list
% loop over channels
delInds = false(1,nBKs);
for k = 1:nBKs
    for m = 1:numel(branch_kinks.channelNums{k})
        chanNum = branch_kinks.channelNums{k}(m);
        delInds(k) = ismemberf(branch_kinks.nurbsParams{k}(m),...
                              chanItrsectParams{chanNum},...
                              'tol',tolParams);
        if delInds(k)
            break
        end
    end
end

% update number of branch points and kinks
nBKs = nnz(~delInds); 
if nBKs < 1
    return
end

% delete branch points or kinks that have been found as intersections
branch_kinks.XX(:,delInds) = [];
branch_kinks.channelNums(delInds) = [];
branch_kinks.ctrlPtNums(delInds) = [];
branch_kinks.nurbsParams(delInds) = [];
new2oldBKnum = find(~delInds);

newNode_coords = nan(nBKs,2);

count = 0;

if isempty(DT)
    error('DT unavailable')
else
    [elemNums,baryCoords] = pointLocation(DT,branch_kinks.XX'); 
end

for k = 1:nBKs
    el = elemNums(k);
    if isnan(el)
        warning('branching point or kink %i cannot be found in any element',k);
        continue
    end
    locEdge = loc_edge_from_bary_coords(baryCoords(k,:),tolBary);
    % Just need to make sure that the kink is included in one of the faces
    % as the info will be distributed to all elements sharing the kink
    % by parent_element_nurbs.m
    if locEdge > 0
        edge = elem_edges(el,locEdge);                    
        fprintf('branch or kink %i on edge %i of elem %i but has not been found\n', ...
                new2oldBKnum(k),edge,el);
        itrsect(edge).num = itrsect(edge).num+1;
        count = count+1;
        newNode_coords(count,:) = branch_kinks.XX(:,k)';
        itrsect(edge).node(end+1) = nNodes+count;
        itrsect(edge).nSeg(end+1) = numel(branch_kinks.channelNums{k});
        itrsect(edge).seg{end+1} = branch_kinks.channelNums{k};
        itrsect(edge).x(end+1) =  branch_kinks.XX(1,k);
        itrsect(edge).y(end+1) =  branch_kinks.XX(2,k);
        itrsect(edge).nurbsParam{end+1} = branch_kinks.nurbsParams{k};
        if calcItrsectVel
            itrsect(edge).vx{end+1} = [1,0];
            itrsect(edge).vy{end+1} = [0,1];
            itrsect(edge).designParamNum{end+1} ...
                = chanDsgnParamNums{branch_kinks.channelNums{k}(1)}(1:2,...
                                          branch_kinks.ctrlPtNums{k}(1))';
        end
    else
         fprintf('branch or kink %i in the interior of elem %i\n',new2oldBKnum(k),el);
         elem_branch_kinks(el).bkNums(end+1) = new2oldBKnum(k);
         count = count+1;
         newNode_coords(count,:) = branch_kinks.XX(:,k)';
         elem_branch_kinks(el).nodes(end+1) = nNodes+count;
    end
end

nodeCoords = [nodeCoords;newNode_coords(1:count,:)];
nNodes = size(nodeCoords,1);
end



