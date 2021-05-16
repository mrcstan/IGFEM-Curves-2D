%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Marcus Tan on 2/2/2015
%%% Copyright 2015 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this function calculates the derivative of each channel length wrt
% design parameters, which can be control points or diameters.
% the derivative wrt the later vanishes
% It also calculates the derivative of each diameter wrt design parameters/
% the derivative wrt control points vanishes 
% ASSUMPTION: Currently, this function is only valid for channels defined 
%             by straight line segments
% OUTPUT:
%   D2L = [DL1/Dd1,DL2/Dd1,...;
%          DL1/Dd2,DL2/Dd2,...;
%                          ...];
%   D2diam = [Ddiam1/Dd1,Ddiam2/Dd1,...;
%             Ddiam1/Dd2,Ddiam2/Dd2,...;
%                          ...];
function [D2L, D2diam] = D2_channel_lengths_n_diameters(channels,designParams)
nChannels = numel(channels.nurbs);
% should form sparse matrix from triplets instead of a full matrix
% but since matrix size is quite small, this shouldn't matter much
D2L = zeros(nChannels,designParams.nParams);
D2diam = zeros(nChannels,designParams.nParams);
for i = 1:designParams.nParams
    if (strcmpi(designParams.type{i},'CTRL_PT'))
        chanNums = designParams.channelNum{i};
        for j = 1:numel(chanNums)
            if(channels.nurbs(chanNums(j)).order == 2)
                ctrlPtDim = designParams.ctrlPtDim(i);
                ctrlPtNum = designParams.ctrlPtNum{i}(j);
                if (ctrlPtNum == 1)
                    %lineSeg = reshape(channels.lineSegs{chanNums(j)}(1,:),[2,2]);
                    %segLength = norm(lineSeg(:,2)-lineSeg(:,1));
                    D2L(chanNums(j),i) = (channels.nurbs(chanNums(j)).coefs(ctrlPtDim,1) ...
                                          -channels.nurbs(chanNums(j)).coefs(ctrlPtDim,2)) ...
                                          /channels.segLengths{chanNums(j)}(1);
                elseif (ctrlPtNum == channels.nurbs(chanNums(j)).number)
                    %lineSeg = reshape(channels.lineSegs{chanNums(j)}(end,:),[2,2]);
                    %segLength = norm(lineSeg(:,2)-lineSeg(:,1));
                    D2L(chanNums(j),i) = (channels.nurbs(chanNums(j)).coefs(ctrlPtDim,ctrlPtNum) ...
                                          -channels.nurbs(chanNums(j)).coefs(ctrlPtDim,ctrlPtNum-1)) ...
                                          /channels.segLengths{chanNums(j)}(end);
                else
                    %lineSeg = reshape(channels.lineSegs{chanNums(j)}(ctrlPtNum-1,:),[2,2]);
                    %segLength1 = norm(lineSeg(:,2)-lineSeg(:,1));
                    %lineSeg = reshape(channels.lineSegs{chanNums(j)}(ctrlPtNum,:),[2,2]);
                    %segLength2 = norm(lineSeg(:,2)-lineSeg(:,1));
                    D2L(chanNums(j),i) = (channels.nurbs(chanNums(j)).coefs(ctrlPtDim,ctrlPtNum) ...
                                          -channels.nurbs(chanNums(j)).coefs(ctrlPtDim,ctrlPtNum-1)) ...
                                          /channels.segLengths{chanNums(j)}(ctrlPtNum-1) ...
                                         +(channels.nurbs(chanNums(j)).coefs(ctrlPtDim,ctrlPtNum) ...
                                           -channels.nurbs(chanNums(j)).coefs(ctrlPtDim,ctrlPtNum+1)) ...
                                          /channels.segLengths{chanNums(j)}(ctrlPtNum);
                                          
                end
                
            else
                error('not yet implemented for general NURBS')
            end
        end
    elseif (strcmpi(designParams.type{i},'DIAM'))
        if (numel(designParams.channelNum{i}) > 1)
            error('only one channel is allowed for diameter type design parameter')
        end
        D2diam(designParams.channelNum{i},i) = 1.0;
    else
        error('unknown design parameters')
    end
end
end
