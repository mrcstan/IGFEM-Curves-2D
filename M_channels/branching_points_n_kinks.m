%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Marcus Tan on 3/7/2016
%%% Copyright 2016 University of Illinois at Urbana-Champaign. 
%%% All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this function finds the branching points and kinks of a channel network
function branch_kinks = branching_points_n_kinks(channels)
branch_kinks = struct('isBranch',[],'XX',[],'channelNums',[], ...
                      'ctrlPtNums',[],'nurbsParams',[]);

% get branching points
[~,ia] = unique(channels.contvty);
duplicatedInds = setdiff(1:numel(channels.contvty),ia);
branchNodes = unique(channels.contvty(duplicatedInds));
branch_kinks.XX = [branch_kinks.XX,channels.pts(branchNodes,:)'];
nBranchNodes = numel(branchNodes);
for i = 1:nBranchNodes
    Lia = ismember(channels.contvty,branchNodes(i));
    channelNums = find(any(Lia,2));
    branch_kinks.channelNums{end+1} = channelNums(:)';
    nChans = numel(channelNums);
    branch_kinks.ctrlPtNums{end+1} = nan(1,nChans);
    branch_kinks.nurbsParams{end+1} = nan(1,nChans);
    for j = 1:nChans
        if Lia(channelNums(j),1)
            branch_kinks.ctrlPtNums{end}(j) = 1;
            branch_kinks.nurbsParams{end}(j) = 0.0;
        else
            branch_kinks.ctrlPtNums{end}(j) ...
                = channels.nurbs(channelNums(j)).number;
            branch_kinks.nurbsParams{end}(j) = 1.0;
        end
    end
end

% get kinks
for i = 1:channels.nNurbs
    if (channels.nurbs(i).order == 2 && channels.nurbs(i).number > 2)
        branch_kinks.XX = [branch_kinks.XX, ...
                           bsxfun(@rdivide, ...
                                  channels.nurbs(i).coefs(1:2,2:end-1), ...
                                  channels.nurbs(i).coefs(4,2:end-1))];                       
        branch_kinks.channelNums ...
            = [branch_kinks.channelNums, ...
                    num2cell(i*ones(1,channels.nurbs(i).number-2))]; 
        branch_kinks.ctrlPtNums ...
            = [branch_kinks.ctrlPtNums, ...
                    num2cell(2:(channels.nurbs(i).number-1))];
        branch_kinks.nurbsParams ...
            = [branch_kinks.nurbsParams, ...
                    num2cell(channels.nurbs(i).knots((channels.nurbs(i).order+1):channels.nurbs(i).number))];
    end
end
branch_kinks.isBranch = false(1,size(branch_kinks.XX,2));
branch_kinks.isBranch(1:nBranchNodes) = true;
end