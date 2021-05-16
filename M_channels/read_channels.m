%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Marcus Tan on 10/11/2014
%%% Last modified date: 10/14/2015
%%% Copyright 2014 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read channel file input
% Note: a keyword and a number is separated by a comma 
%       Other than the newline character or a blank, no other separator is
%       used
% INPUT: 
%  inputFile:
%  options.figBounds: for plotting the bounding boxes of the control point design parameters
%  options.filePrefix: the prefix of the filenames of the output figures
%  options.guimode:
%  options.boundsFile:
%  options.maxTrials: max number of trials to generate random configuration
%                     that satisfies nonlinear constraints
%  options.nlconfun: a function handle to the nonlinear constraint
%  options.nlcon: other inputs for nonlinear constraint including the rhs
%  options.generatePolygons: generate polygons that are needed by some
%                            nonlinear constraints
% OUTPUT:
%  channels:
%  designParams:
% FORMAT:
% # A comment line starts with a hex. The rest of the line will then be
% # ignored
% number of channels, <N>
%
% end point number and temperature, <npres>
% i_1, T_1
% ...
% i_npres, T_npres
%
% model, <model type>
%
% EITHER  
% mcf, <nch>
% mcf_1,
% ...
% mcf_nch
% OR
% mass in, <mass flow rate in>
%
% heat capacity, <heat capacity>
%
% viscosity, <viscosity>
%
% pressure out, <pressure out>
%
% diameters, <nch>
%
% nurbs control points
% xc_11 ... xc_nd1 w1
% xc_12 ... xc_nd2 w2
% ...
% xc_1ncpts ... xc_ndncpts wncpts
% nurbs knot vector
% xi1 ... xi_knots
% connectivity, 1
% i j 
% <empty line>
% nurbs control points
% xc_11 ... xc_nd1 w1
% xc_12 ... xc_nd2 w2
% ...
% xc_1ncpts ... xc_ndncpts wncpts
% nurbs knot vector
% xi1 ... xi_knots
% connectivity
% i j 
% 
% number of design variables, <number of design variables>
% design
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
% nd: number of dimensions (nd <= 2)
% npts: number of channel end points
% nch: number of channels
% npres: number of channel end points where temperature is prescribed
% ncpts: number of control points
% model type: 'mean temperature' or 'constant heat flux'
% REMARKS:
%   i) If nd == 4, the control points are assumed not yet multiplied by the
%       weights
%--------------------------------------------------------------------------
% EXAMPLE:
%{
number of channels, 7

end point number and temperature, 1
1 21.5

model, mean temperature
mass in, 1, 5e-4
heat capacity, 3494
viscosity, 3.405e-6
pressure out, 8, 0
cross section, square 

diameters
7.5e-4
7.5e-4
7.5e-4
7.5e-4
7.5e-4
7.5e-4
7.5e-4

nurbs control points
0 0.18
0.02 0.18
nurbs knot vector
0 0 1 1
connectivity, 1
1 2

nurbs control points
0.02 0.18
0.02 0.1
nurbs knot vector
0 0 1 1
connectivity
2 4

nurbs control points
0.02 0.18
0.13 0.18
0.13 0.1
nurbs knot vector
0 0 0.5 1 1
connectivity
2 5

nurbs control points
0.02 0.1
0.13 0.1
nurbs knot vector
0 0 1 1
connectivity
4 5

nurbs control points
0.02 0.1
0.02 0.02
0.13 0.02
nurbs knot vector
0 0 0.5 1 1
connectivity
4 7

nurbs control points
0.13 0.1
0.13 0.02
nurbs knot vector
0 0 1 1
connectivity
5 7

nurbs control points
0.13 0.02
0.15 0.02
nurbs knot vector
0 0 1 1
connectivity
7 8

number of design variables, 12
design
channel, 1, cpt, 2, x, 0.005, 0.07
design
channel, 1, cpt, 2, y, 0.145, 0.195
design
channel, 3, cpt, 2, x, 0.08, 0.145
design
channel, 3, cpt, 2, y, 0.145, 0.195
design
channel, 4, cpt, 1, x, 0.005, 0.07
design
channel, 4, cpt, 1, y, 0.065, 0.135
design
channel, 4, cpt, 2, x, 0.08, 0.145
design
channel, 4, cpt, 2, y, 0.065, 0.135
design
channel, 5, cpt, 2, x, 0.005, 0.07
design
channel, 5, cpt, 2, y, 0.005, 0.055
design
channel, 5, cpt, 3, x, 0.08, 0.145
design
channel, 5, cpt, 3, y, 0.005, 0.055  
%}
function [channels,designParams] = read_channels(inputFile)                                   

if (~ischar(inputFile))
    error('input file name must be a character array');
end

% channel configuration
nChannelKeyword = 'number of channels';
mcfKeyword = 'mcf';
endPtNTempKeyword = 'end point number and temperature';
modelKeyword = 'model';
massInKeyword = 'mass in';
powerXdensityKeyword = 'powerXdensity';
heatCapacityKeyword = 'heat capacity';
densityKeyword = 'density';
defaultDensity = 1065;
viscosityKeyword = 'viscosity'; % kinematic viscosity
visModelKeyword = 'viscosity model'; % kinematic viscosity
defaultVisModel = 1; % temperature dependent viscosity (see kinematic_viscosity.m for more info)
pressureOutKeyword = 'pressure out';
crossSectionKeyword = 'cross section';
validCrossSections = {'square','circular','rectangular'};
diameterKeyword = 'diameters';
heightKeyword = 'heights';
nurbsCtrlPtKeyword = 'nurbs control points';
nurbsKnotKeyword = 'nurbs knot vector';
connectivityKeyword = 'connectivity';

% design variables
nDesignKeyword = 'number of design variables';
designKeyword = 'design';
channelKeyword = 'channel';
ctrlPtKeyword = 'cpt';
ctrlPtDimKeyword = {'x','y','z'};
diamKeyword = 'diam';

ndesignParams = 0;

channels = struct;

modelType = [];
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

    if (strcmpi(splitStr{1},nChannelKeyword))       
        nChannels = str2double(splitStr{2});
        if (isnan(nChannels))
            error('number of channels Keyword should be followed by comma, number of channels')
        end
    elseif (strcmpi(splitStr{1},endPtNTempKeyword))
        nPrescribedEndPts = str2double(splitStr{2});
        if (isnan(nPrescribedEndPts))
            error('channel connectivityend point number and temperature Keyword should be followed by comma and number of prescribed end points')
        end
        channels.pt_temp = fscanf(fileID,'%f %f',[2,nPrescribedEndPts])';
    elseif (strcmpi(splitStr{1},modelKeyword))
        modelType = strtrim(splitStr{2});
    elseif (strcmpi(splitStr{1},massInKeyword))
        if isfield(channels,'inletEndPoint')
            % handles case of multiple inlets
            channels.inletEndPoint(end+1) =  str2double(splitStr{2});
        else
            channels.inletEndPoint = str2double(splitStr{2});
        end
        if isfield(channels,'massin')
            % handles case of multiple inlets
            channels.massin(end+1) = str2double(splitStr{3});
        else
            channels.massin = str2double(splitStr{3});
        end 
    elseif (strcmpi(splitStr{1},powerXdensityKeyword))
        channels.inletEndPoint = str2double(splitStr{2});
        channels.powerXdensity = str2double(splitStr{3});
    elseif (strcmpi(splitStr{1},heatCapacityKeyword))
        channels.heatCapacity = str2double(splitStr{2});
    elseif (strcmpi(splitStr{1},densityKeyword))
        channels.density = str2double(splitStr{2});
    elseif (strcmpi(splitStr{1},viscosityKeyword))
        channels.viscosity = str2double(splitStr{2});
    elseif (strcmpi(splitStr{1},visModelKeyword))
        channels.viscosityModel = str2double(splitStr{2});
    elseif (strcmpi(splitStr{1},pressureOutKeyword))
        if isfield(channels,'pressureOutletEndPoint')
             % handles case of multiple outlets
            channels.pressureOutletEndPoint(end+1) = str2double(splitStr{2});
        else
            channels.pressureOutletEndPoint = str2double(splitStr{2});
        end
        if isfield(channels,'pressureOut')
             % handles case of multiple outlets
            channels.pressureOut(end+1) = str2double(splitStr{3});
        else
            channels.pressureOut = str2double(splitStr{3});
        end
    elseif (strcmpi(splitStr{1},crossSectionKeyword))
        splitStr{2} = strtrim(splitStr{2});
        if any(strcmpi(splitStr{2},validCrossSections))
            channels.crossSection = splitStr{2};
        else
            error('invalid cross section')
        end
    elseif (strcmpi(splitStr{1},diameterKeyword))
        nChannelDiams = 0;
        channels.diams = [];
        line = fgetl(fileID);
        num = str2double(line);
        while ~isnan(num)
            nChannelDiams = nChannelDiams+1;
            channels.diams = [channels.diams;num];
            line = fgetl(fileID);
            num = str2double(line);
        end
    elseif (strcmpi(splitStr{1},heightKeyword))
        nChannelHeights = 0;
        channels.heights = [];
        line = fgetl(fileID);
        num = str2double(line);
        while ~isnan(num)
            nChannelHeights = nChannelHeights+1;
            channels.heights = [channels.heights;num];
            line = fgetl(fileID);
            num = str2double(line);
        end
    elseif (strcmpi(splitStr{1},mcfKeyword))
        nChannelsMcf = 0;
        channels.mcf = [];
        line = fgetl(fileID);
        num = str2double(line);
        while ~isnan(num)
            nChannelsMcf = nChannelsMcf + 1;
            channels.mcf = [channels.mcf;num];
            line = fgetl(fileID);
            num = str2double(line);
        end
    elseif (strcmpi(splitStr{1},nDesignKeyword))
        ndesignParams = str2double(splitStr{2});
        if (isnan(ndesignParams))
            error('number of design variables keyword should be followed by comma and the number of design variables')
        end
    end
    line = fgetl(fileID);
end

if (isfield(channels,'mcf') && nChannels ~= nChannelsMcf)
    error('number of provided mcf is different from the number of channels')
end

if ~isfield(channels,'massin') && ~isfield(channels,'powerXdensity')
    error('either mass in or powerXdensity need to be specified')
elseif isfield(channels,'massin') && isfield(channels,'powerXdensity')
    error('specify either mass in or powerXdensity')
elseif isfield(channels,'massin')
    channels.powerXdensity = [];
else
    channels.massin = [];
end
    
    
if (isfield(channels,'diams') && nChannels ~= nChannelDiams) 
    error('number of provided channel diameters is different from the number of channels')
end

if (isfield(channels,'heights') && nChannels ~= nChannelHeights) 
    error('number of provided channel heights is different from the number of channels')
end


if (~isfield(channels,'crossSection'))
    warning('cross section of channel not specified, use %s as default',...
            validCrossSections{1})
    channels.crossSection = validCrossSections{1};
end

if ~isfield(channels,'heights')
    if strcmpi(channels.crossSection,validCrossSections{3})
        error('channel heights must be provided for rectangular cross sections')
    else
        channels.heights = [];
    end
end

if (~isfield(channels,'density'))
    warning('density of fluid not specified, use %g as default ', defaultDensity)
    channels.density = defaultDensity;
end

if (~isfield(channels,'viscosityModel'))
    warning('viscosity model not specified, use model %i as default ', defaultVisModel)
    channels.viscosityModel = defaultVisModel;
end

if (ndesignParams == 0)
    warning('no design parameters specified')
end

frewind(fileID); % Back to the beginning of the input file
line = fgetl(fileID);
p = cell(nChannels);
knots= cell(nChannels);
channels.contvty = nan(nChannels,2);

designParams.nParams = ndesignParams;
designParams.channelNum = cell(ndesignParams,1);
designParams.ctrlPtNum = cell(ndesignParams,1);
designParams.ctrlPtDim = nan(ndesignParams,1);
designParams.type = cell(ndesignParams,1);
designParams.bounds = inf(ndesignParams,2);

channelNum = 0;
designNum = 0;
while ~feof(fileID)
    line = strtrim(line);
     if (isempty(line) ||strncmpi(line,'#',1))
        line = fgetl(fileID);
        continue
    end
    splitStr = regexp(line,',','split');
    splitStr{1} = strtrim(splitStr{1});
    while strcmpi(splitStr{1},nurbsCtrlPtKeyword)
        channelNum = channelNum + 1;  
        line = strtrim(fgetl(fileID));
        nums = str2num(line);
        while ~isempty(nums)
            p{channelNum} = [p{channelNum},nums'];
            line = strtrim(fgetl(fileID));
            nums = str2num(line);
        end
        
        splitStr = regexp(line,',','split');
        splitStr{1} = strtrim(splitStr{1});
        % when another nurbsCtrlPtKeyword is met, the subsequent lines
        % belong to a new channel
        while ~strcmpi(splitStr{1},nurbsCtrlPtKeyword)
            line = strtrim(line);
            splitStr = regexp(line,',','split');
            splitStr{1} = strtrim(splitStr{1});
            if (strcmpi(splitStr{1},nurbsKnotKeyword))                
                line = strtrim(fgetl(fileID));
                knots{channelNum} = str2num(line);
            elseif (strcmpi(splitStr{1},connectivityKeyword))
                line = strtrim(fgetl(fileID));
                channels.contvty(channelNum,:) = str2num(line);
            end            
            if feof(fileID)
                break
            end
            line = fgetl(fileID);
            splitStr = regexp(line,',','split');
            splitStr{1} = strtrim(splitStr{1});
            if (isempty(splitStr{1}) ...
                || strcmpi(splitStr{1},designKeyword))
                break
            end
        end
    end
 
    if (strcmpi(splitStr{1},designKeyword))
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

if designNum ~= ndesignParams
    error('number of provided design parameters not equal to specified number of design parameters')
end
channels.nNurbs = nChannels; % number of nurbs

% compact the channel end point numbers so that they are numbered consecutively
% this result in a new numbering for the channel end points
[oldNodes,~,newContvty] = unique(channels.contvty(:));
% determine if the nodes have been compacted
if (any(oldNodes' ~= 1:numel(oldNodes)))
    warning('channel node numbers have been compacted') 
    compactNodes = true;
     channels.contvty = reshape(newContvty,numel(newContvty)/2,2);
else
    compactNodes = false;
end
channels.pts = nan(max(channels.contvty(:)),2);

for i = 1:nChannels
    if (isempty(p{i}) || isempty(knots{i}) || any(isnan(channels.contvty(i,:)))) 
        error('the control points, the knots or the connectivity a channel is not specified')
    end
    if size(p{i},1) == 4 && (numel(knots{i})-size(p{i},2)) == 2
        error('linear NURBS with non-unity weights not supported')
    end
    if size(p{i},1) == 4
        p{i}(1,:) = p{i}(1,:).*p{i}(4,:);
        p{i}(2,:) = p{i}(2,:).*p{i}(4,:);
        p{i}(3,:) = p{i}(3,:).*p{i}(4,:);
    end
    channels.nurbs(i) = nrbmak(p{i},knots{i});
    channels.pts(channels.contvty(i,1),:) = [p{i}(1,1),p{i}(2,1)];
    channels.pts(channels.contvty(i,2),:) = [p{i}(1,end),p{i}(2,end)];
end

if (isfield(channels,'massin') && isfield(channels,'mcf'))
    warning('mass in takes precedence when both mcf and mass in are specified')
end
if compactNodes
    warning('inlets and outlets have been relabeled')
    % change inletEndPoint to new the new numbering of the channel
    % end points
    for i = 1:numel(channels.inletEndPoint)
        channels.inletEndPoint(i) ...
            = find(oldNodes == channels.inletEndPoint(i), 1,'first');
    end
    for i = 1:numel(channels.pressureOutletEndPoint)
        channels.pressureOutletEndPoint(i) ...
            = find(oldNodes == channels.pressureOutletEndPoint,1,'first');
    end
end 

% modelType: 1 - dT/ds model
%            2 - h(T(s)-Tin) model
if (strcmpi(modelType,'mean temperature'))
    channels.model = 1;
elseif (strcmpi(modelType,'constant heat'))
    channels.model = 2; 
else
    error('unrecognized model type')
end


% for mex function to find intersections
channels.kind = ones(channels.nNurbs,1);
for i = 1:channels.nNurbs
    if (any(channels.nurbs(i).coefs(4,:) ~= 1))
        channels.kind(i) = 2;
        break
    end
end
channels.rows{1} = [1,2]; % 2D B-splines
channels.rows{2} = [1,2,4]; % 2D NURBS

%channels = rmfield(channels,'length');


if (channels.model == 2)
    channels.pt_temp = [];
end

end
