%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Marcus Tan on 10/11/2014
%%% Last modified date: 1/27/2016
%%% Copyright 2014 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read channel file input that contains bounds on design parameters and
% randomize the control points within the bounds and accept only those that 
% satisfy a nonlinear constraint nlconfun
% Note: a keyword and a number is separated by a comma 
%       Other than the newline character or a blank, no other separator is
%       used
% INPUT: (if design parameters are not specified in input file only first
%         variables is relevant)
%  inputFile:
% FORMAT:
% number of design variables, <number of design variables>
% design, <design variable number>
% channel, <channel number>, <design variable type>, other keywords and
% parameters depending on <design variable type>
% if <design variable type> == 'cpt', 'cpt,' must be followed by
% <control point number>, <x or y or z>, <lower bound>, <upper bound>
% else if <design variable type> == 'diam', <lower bound>, <upper bound>
% else
% end
% set <lower bound> to -inf or <upper bound> to inf if unbounded below or
% above
% 
% FORMAT variable names and possible values:
% nd: number of dimensions
% npts: number of channel end points
% nch: number of channels
% npres: number of channel end points where temperature is prescribed
% ncpts: number of control points
% model type: 'mean temperature' or 'constant heat flux'
%--------------------------------------------------------------------------
% EXAMPLE:
%{
number of design variables, 12
design, 1
channel, 1, cpt, 2, x, 0.005, 0.07
design, 2
channel, 1, cpt, 2, y, 0.145, 0.195
design, 3 
channel, 3, cpt, 2, x, 0.08, 0.145
design, 4 
channel, 3, cpt, 2, y, 0.145, 0.195
design, 5
channel, 4, cpt, 1, x, 0.005, 0.07
design, 6
channel, 4, cpt, 1, y, 0.065, 0.135
design, 7
channel, 4, cpt, 2, x, 0.08, 0.145
design, 8
channel, 4, cpt, 2, y, 0.065, 0.135
design, 9
channel, 5, cpt, 2, x, 0.005, 0.07
design, 10
channel, 5, cpt, 2, y, 0.005, 0.055
design, 11
channel, 5, cpt, 3, x, 0.08, 0.145
design, 12
channel, 5, cpt, 3, y, 0.005, 0.055  
%}
function [channels,designParamIniVals,restrictedParamIniVals] ...
                = randomize_channels(channels,...
                                     vertices2params,...
                                     restrictedParams,...
                                     inputFile,...
                                     maxTrials,...
                                     nlconfun,...
                                     nlcon)

if (~ischar(inputFile))
    error('input file name must be a character array');
end

% design variables
nDesignKeyword = 'number of design variables';
designKeyword = 'design';
channelKeyword = 'channel';
ctrlPtKeyword = 'cpt';
ctrlPtDimKeyword = {'x','y','z'};
diamKeyword = 'diam';
ndesignParams = 0;

fileID = fopen(inputFile,'r');
line = fgetl(fileID);

while (ischar(line))
    line = strtrim(line);  
    % skip comment line (starting with #) and empty line
    if (isempty(line) ||strncmpi(line,'#',1))
        line = fgetl(fileID);
        continue
    end
    splitStr = regexp(line,',','split');
    splitStr{1} = strtrim(splitStr{1});
    if (strcmpi(splitStr{1},nDesignKeyword))
        ndesignParams = str2double(splitStr{2});
        if (isnan(ndesignParams))
            error('number of design variables keyword should be followed by comma and the number of design variables')
        end
    end
    line = fgetl(fileID);
end


if (ndesignParams == 0)
    error('no design parameters specified')
end

frewind(fileID); % Back to the beginning of the input file
line = fgetl(fileID);

designParams.nParams = ndesignParams;
designParams.channelNum = cell(ndesignParams,1);
designParams.ctrlPtNum = cell(ndesignParams,1);
designParams.ctrlPtDim = nan(ndesignParams,1);
designParams.type = cell(ndesignParams,1);
designParams.bounds = inf(ndesignParams,2);
designParams.vertices2params = vertices2params;

designNum = 0;
while ~feof(fileID)
    line = strtrim(line);
    if (isempty(line) ||strncmpi(line,'#',1))
        line = fgetl(fileID);
        continue
    end
    splitStr = regexp(line,',','split');
    splitStr{1} = strtrim(splitStr{1});
     if strcmpi(splitStr{1},designKeyword)
        designNum = designNum + 1;
        line = fgetl(fileID);
        splitStr = regexp(line,',','split');
        for i = 1:numel(splitStr)
            splitStr{i} = strtrim(splitStr{i});
        end
        if (strcmpi(splitStr{1},channelKeyword))
            designParams.channelNum{designNum} = str2double(splitStr{2});  
            %chanNum = designParams.channelNum{designNum};
            if (strcmpi(splitStr{3},ctrlPtKeyword))
                designParams.ctrlPtNum{designNum} = str2double(splitStr{4});
                if (isnan(designParams.ctrlPtNum{designNum}))
                    error('control point keyword should be followed by comma and control point number')
                end
                if (strcmpi(splitStr{5},ctrlPtDimKeyword{1}))
                    designParams.ctrlPtDim(designNum) = 1;
                elseif (strcmpi(splitStr{5},ctrlPtDimKeyword{2}))
                    designParams.ctrlPtDim(designNum) = 2;
                elseif (strcmpi(splitStr{5},ctrlPtDimKeyword{3}))
                    designParams.ctrlPtDim(designNum) = 3;
                else
                    error('unrecognized dimension of control point')
                end
                
                %dim = designParams.ctrlPtDim(designNum);
                %num = designParams.ctrlPtNum{designNum};
                %designParams.iniVals(designNum) ...
                %    = channels.nurbs(chanNum).coefs(dim,num);
                designParams.bounds(designNum,1) = str2double(splitStr{6});
                designParams.bounds(designNum,2) = str2double(splitStr{7});
                designParams.type{designNum} = 'CTRL_PT';
                
            elseif (strcmpi(splitStr{3},diamKeyword))
                %designParams.iniVals(designNum) = channels.diam(chanNum);
                designParams.bounds(designNum,1) = str2double(splitStr{4});
                designParams.bounds(designNum,2) = str2double(splitStr{5});
                designParams.type{designNum} = 'DIAM';
            else
                error('unrecognized design variable type')
            end
        else
            error('channel keyword following design %i not found', designNum)
        end
    end
    
    line = fgetl(fileID);
end
fclose(fileID);

[designParams,channels.designParamNum] ...
    = design_params2channel_params_map(designParams,...
                                       channels);

generated = false;
% randomize the initial values                
try
    del = zeros(designParams.nParams+restrictedParams.nParams,1);
    for i = 1:maxTrials
        designParams.iniVals = randomized_bounded_values(designParams.bounds);
        restrictedParams.iniVals ...
            = update_restricted_params_ini_vals(designParams, ...
                                                restrictedParams);
        channels = update_channels([designParams.iniVals;restrictedParams.iniVals],...
                                   designParams, ...
                                   restrictedParams, ...
                                   channels, ...
                                   'replace', ...
                                   false);
        if (any(nlconfun(del,designParams,restrictedParams,channels,nlcon) > 0))
            continue
        end
        fprintf('number of attempts to get random channel = %i \n',i)
        generated = true;
        break;
    end
   
catch err
    warning(getReport(err));
  
    designParams.iniVals= randomized_bounded_values(designParams.bounds);
    channels = update_channels(designParams.iniVals, ...
                               designParams, ...
                               restrictedParams,...
                               channels, ...
                               'replace', ...
                               false);
    restrictedParams.iniVals ...
            = update_restricted_params_ini_vals(designParams, ...
                                                restrictedParams);                       
    generated = true;
    
end

if (~generated)
    error('maximum number of attempts to randomize channel has been exceeded')
else
    designParamIniVals = designParams.iniVals;
    restrictedParamIniVals = restrictedParams.iniVals;
end                       

end
