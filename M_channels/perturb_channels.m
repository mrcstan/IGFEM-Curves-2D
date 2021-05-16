%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Marcus Tan on 12/21/2014
%%% Copyright 2014 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this function perturbs a control point of a channel closest to coord  
% by a given amount and
%  update the other channels if they share the same control point, i.e., if
% the control point is a branching point
% ASSUMPTION: All weights are unity
function channels = perturb_channels(channelNum, channels, coord, del)
% distance of each control point from coord
distsq = (channels.nurbs(channelNum).coefs(1,:)-coord(1)).^2 ...
        +(channels.nurbs(channelNum).coefs(2,:)-coord(2)).^2; 
[~,indMin] = min(distsq); % control point number of control point closest to coord
channels.nurbs(channelNum).coefs(1:2,indMin) ...
    = channels.nurbs(channelNum).coefs(1:2,indMin) + del';
if (indMin ~= 1  && indMin ~= channels.nurbs(channelNum).number)
    return
end

% if perturbed control point is an end point of the channel, channels
% connected to this control point are affected too and they must be updated
% accordingly
if (indMin == 1)
    endPtInd = 1;
else
    endPtInd = 2;
end
ctrlPtCoord = channels.nurbs(channelNum).coefs(1:2,indMin);
channels.pts(channels.contvty(channelNum,endPtInd),:) = ctrlPtCoord';
for i = 1:size(channels.contvty,1)        
    if (i ~= channelNum)
        endPtInd2 = find(channels.contvty(i,:) == channels.contvty(channelNum,endPtInd));
        if (endPtInd2 == 1)
            channels.nurbs(i).coefs(1:2,1) = ctrlPtCoord;
        elseif (endPtInd2 == 2)
            channels.nurbs(i).coefs(1:2,end) = ctrlPtCoord;
        end
    end
end

end