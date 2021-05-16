%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Marcus Tan on 11/4/2014
%%% Last modified date: 12/21/2014
%%% Copyright 2014 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this function maps the design parameter numbers to the channel parameters
% OUTPUT:
%   channelDesignParamNum: an array of cells of length = number of channels
%                          each cell corresponds to a channel and contains
%                          a 3xnumber of control points matrix.
%                          each entry contains the design parameter number
%                          of the corresponding coordinate if it is chosen
%                          as a design parameter. otherwise, the entry is
%                          nan. the first, second and third rows correspond
%                          to x-, y-coordinates and weights of the control points
%   designParams:
%       designParams.iniVals:
%       designParams.channelNum:
%       designParams.ctrlPtNum:
function [designParams,channelDesignParamNum] ...
                = design_params2channel_params_map(inDesignParams,...
                                                   channels)
if (nargout > 1)
    nChannels = numel(channels.nurbs);
    channelDesignParamNum = cell(nChannels,1);
    for i = 1:nChannels
        channelDesignParamNum{i} = nan(3,channels.nurbs(i).number); 
    end
end
designParams = inDesignParams;
designParams.iniVals = nan(designParams.nParams,1);
for i = 1:inDesignParams.nParams
    channelNum = inDesignParams.channelNum{i};
    % check validity of design parameter channel number
    if (channelNum < 1 || channelNum > numel(channels.nurbs))
        error('design parameter %i channel number out of bound',i)
    end
    if (strcmpi(inDesignParams.type{i},'CTRL_PT'))
        % check validity of control point design parameter number
        if (inDesignParams.ctrlPtDim(i) < 1 || inDesignParams.ctrlPtDim(i) > 4)
            error('control point design parameter %i dimension', i) 
        end
        if (inDesignParams.ctrlPtNum{i} < 1  ...
                || inDesignParams.ctrlPtNum{i} > channels.nurbs(channelNum).number)
            error('design parameter %i control point number out of bound', i) 
        end
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
            if (nargout > 1)
                channelDesignParamNum{channelNum}(designParams.ctrlPtDim(i),...
                                                  designParams.ctrlPtNum{i}) ...
                                                  = i;
            end
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
                if (nargout > 1)
                    channelDesignParamNum{sharedChannels(j)}(designParams.ctrlPtDim(i),...
                                                             designParams.ctrlPtNum{i}(j)) ...
                                                            = i;
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