%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Marcus Tan: 10/20/2015
%%% Copyright 2015 University of Illinois at Urbana-Champaign. 
%%% All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Note: Tested on MATLAB 2015b
%       graph and maxflow functions are not available in MATLAB 2014a or older
% This function calculates the average channel network redundancy based on
% edge-disjoint paths
function ave = average_network_redundancy(channels,excludedNodes)
    connectivity = channels.contvty;
    %
    % nNodes = max(connectivity(:));
    delInds = nan(size(connectivity,1),1);
    count = 0;
    % delete channels containing excluded nodes
    for i = 1:numel(excludedNodes)
        [rows,~] = find(connectivity == excludedNodes(i));
        delInds(count+1:count+numel(rows)) = rows;
        count = count+numel(rows);
    end
    delInds(isnan(delInds)) = [];
    connectivity(delInds,:) = [];
    %
    nUsedChannels = size(connectivity,1);
    %G = sparse([connectivity(:,1);connectivity(:,2)],...
    %           [connectivity(:,2);connectivity(:,1)],...
    %           ones(2*nUsedChannels,1),nNodes,nNodes);
    G = graph(connectivity(:,1),connectivity(:,2),ones(nUsedChannels,1));
    ave = 0;
    usedNodes = unique(connectivity(:)); 
    usedNodePairs = combnk(usedNodes,2);
    nNodePairs = size(usedNodePairs,1);
    for i = 1:nNodePairs
        ave = ave + maxflow(G,usedNodePairs(i,1),usedNodePairs(i,2));
    end
    ave = ave/nNodePairs;
end