%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Marcus Tan on 11/4/2014
%%% Last modified date: 11/4/2014
%%% Copyright 2014 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this function maps the design parameter numbers to the channel parameters
function designParams = design_params2channel_params_map(inDesignParams,...
                                                         channels)
designParams = inDesignParams;
designParams.iniVals = nan(designParams.nParams,1);
for i = 1:inDesignParams.nParams
    channelNum = inDesignParams.channelNum{i};
    if (strcmpi(inDesignParams.type{i},'CTRL_PT'))
        % if control points are the terminal ends of a nurbs channel,
        % get the numbers as well as the control point numbers of all the
        % channel sharing the point

        designParams.iniVals(i) ...
            = channels.nurbs(channelNum).coefs(inDesignParams.ctrlPtDim(i),...
                                               inDesignParams.ctrlPtNum{i});
        if (inDesignParams.ctrlPtNum{i} == 1)
            endPtInd = 1;
        elseif (inDesignParams.ctrlPtNum{i} == channels.nurbs(channelNum).number)
            endPtInd = 2;
        else
            endPtInd = [];
        end
        if (~isempty(endPtInd))
            [sharedChannels,cols] ...
                = find(channels.contvty == channels.contvty(channelNum,endPtInd));
            nSharedChannels = numel(sharedChannels); 
            designParams.channelNum{i} = zeros(nSharedChannels,1);
            designParams.ctrlPtNum{i} = ones(nSharedChannels,1);
            %designParams.ctrlPtDim{i} = repmat(inDesignParams.ctrlPtDim{i},...
            %                                  [nSharedChannels,1]);
                               
            for j = 1:nSharedChannels
                designParams.channelNum{i}(j) = sharedChannels(j);
                if (cols(j) == 2)
                    designParams.ctrlPtNum{i}(j) ...
                        = channels.nurbs(sharedChannels(j)).number; 
                end
            end
        end
    elseif (strcmpi(inDesignParams.type{i},'DIAM'))
        designParams.iniVals(i) = channels.diams(channelNum);
    else
        error('design parameter number %i is not defined',i)
    end
end
end