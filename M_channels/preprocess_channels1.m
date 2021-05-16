%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Marcus Tan on 10/11/2014
%%% Last modified date: 2/5/2015
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
%  options.boundsFile: file with same format as channel file containing 
%                       the bounds of each control point
%                       design parameters for randomization.
%  options.lhsFile: Latin hypercube sampling file. if both boundsFile and
%                   lhsFile are present, lhsFile is used.
%  options.maxTrials: max number of trials to generate random configuration
%                     that satisfies nonlinear constraints
%  options.nlconfun: a function handle to the nonlinear constraint
%  options.nlcon: other inputs for nonlinear constraint including the rhs
%  options.polygonGenFile: file with same format as channel file to generate 
%                       polygons that are needed by some
%                       nonlinear constraints. this is guaranteed to work 
%                       only for reference configuration
%  options.polygonFileName: polygons output file name
% OUTPUT:
%  channels:
%  bodySource:
%  designParams:
%  restrictedParams:
function varargout = preprocess_channels(inputFile,options)                                   
% boundaries of domain
bodySource = [];

if (~ischar(inputFile))
    error('input file name must be a character array');
end

if (nargin < 2)
    options.figBounds = [];
    options.filePrefix = [];
    options.guimode = false;
    options.boundsFile = [];
    options.lhsFile = [];
    options.polygonGenFile = []; 
    options.maxTrials = 50;
    options.nlconfun = [];
    options.nlcon = [];
    
    options.polygonFileName = [];
end
if (~isfield(options,'figBounds'))
    options.figBounds = [];
end
if (~isfield(options,'filePrefix'))
    options.filePrefix = [];
end
if (~isfield(options,'guimode'))
    options.guimode = false;
end
if (~isfield(options,'boundsFile'))
    options.boundsFile = [];
end
if (~isfield(options,'lhsFile'))
    options.lhsFile = [];
end
if (~isfield(options,'maxTrials'))
    options.maxTrials = 50;
end
if (~isfield(options,'nlconfun'))
    options.nlconfun = [];
end
if (~isfield(options,'nlcon'))
    options.nlcon = [];
end
if (~isfield(options,'polygonGenFile'))
    options.polygonGenFile = [];
end

if (~isfield(options,'polygonFileName'))
    options.polygonFileName = [];
end


if (ischar(options.lhsFile) && ischar(options.boundsFile))
    warning('Latin hypercube sampling file is used instead of the bounds file')
    options.boundsFile = [];
end

if (ischar(options.lhsFile))
    if(~isfield(options,'sampleNum'))
        error('must specify sample number to select from Latin hypercube sampling file')
    end
    if(~isnumeric(options.sampleNum))
        error('sample number must be a numeric')
    end
end

[channels,designParams] = read_channels(inputFile);

    

% design parameters
if (~options.guimode && designParams.nParams)
    [designParams,channels.designParamNum] ...
        = design_params2channel_params_map(designParams,...
                                           channels);
end

if (ischar(options.polygonGenFile))
    warning('channel_polygons guaranteed to work only when applied to reference parallel channels')
    [polyChannels,polyDesignParams] = read_channels(options.polygonGenFile);
    [polyDesignParams,polyChannels.designParamNum] ...
        = design_params2channel_params_map(polyDesignParams,...
                                           polyChannels);
    [channels.polygons, channels.vertexCoords,...
     designParams.vertices2params,...
     restrictedParams.nParams, ...
     restrictedParams.iniVals, ...
     restrictedParams.paramPairs] ...
        = channel_polygons(polyChannels,polyDesignParams.nParams,options);
    % make sure that the restricted params have the same values as the
    % corresponding params. 
    restrictedParams.iniVals ...
        = update_restricted_params_ini_vals(designParams, ...
                                            restrictedParams);
else
    restrictedParams.nParams = 0;
    restrictedParams.iniVals = [];
end
% randomize the initial values 
if (ischar(options.boundsFile))
    [channels,designParams.iniVals,restrictedParams.iniVals] ...
                = randomize_channels(channels,...
                                    designParams.vertices2params,...
                                    restrictedParams,...
                                    options.boundsFile,...
                                    options.maxTrials,...
                                    options.nlconfun,...
                                    options.nlcon);
                                                            
end
if (ischar(options.lhsFile))
    designParams.iniVals = read_lhs_file(options.lhsFile,options.sampleNum);
    restrictedParams.iniVals ...
        = update_restricted_params_ini_vals(designParams, ...
                                            restrictedParams);
end
% need to calculate pressure, mass and mcf even for randomFlag == true
% because the update_channels function does not calculate those when those
% output are not requested
if (~options.guimode)
    if (designParams.nParams)
        [channels,pressure,mass] ...
                = update_channels([designParams.iniVals;restrictedParams.iniVals], ...
                                  designParams, ...
                                  restrictedParams, ...
                                  channels, ...
                                  'replace');
    else
        ngpts = 10;
        [channels.length,channels.lineSegs,channels.segLengths]...
                    = nurbs_channel_lengths(channels,ngpts);    
        channels.vol = channel_volume(channels.length,channels.diams,channels.crossSection);
        channels.junc = branching_points(channels.contvty);   
        channels.kinks = kinks(channels.nurbs);
        channels.designParamNum = [];
        [pressure,mass] = network_pressure_mass_flow_rate(channels.contvty,...
                                                          channels.nurbs,...
                                                          channels.diams,...
                                                          channels.viscosity,...
                                                          channels.inletEndPoint,...
                                                          channels.massin,...
                                                          channels.pressureOutletEndPoint,...
                                                          channels.pressureOut,...
                                                          channels.crossSection);
        channels.mcf = mass*channels.heatCapacity;                                               
    end
       
end


% plot pressure and mass flow rate in terms of kPa and ml/min, respectively
if (~options.guimode)
    if (ischar(options.filePrefix))
        fig = figure('visible','off');
        %set(fig,'ResizeFcn','set(gcf,''visible'',''on'')'); 
        plot_channel_network(channels.contvty,channels.nurbs,pressure/1000.0,mass*1e6*60/1065)
        fname = [options.filePrefix,'_ini_mass_pressure'];
        set(fig, 'PaperPositionMode', 'auto');
        print('-djpeg ', '-r300', fname);
        saveas(fig, fname, 'fig');
    else
        figure
        plot_channel_network(channels.contvty,channels.nurbs,pressure/1000.0,mass*1e6*60/1065)
    end
end
if (~options.guimode && designParams.nParams)
    if (ischar(options.filePrefix))
        fig = figure('visible','off');
        %set(fig,'ResizeFcn','set(gcf,''visible'',''on'')');    
        plot_channels_and_design_parameters(channels.contvty,channels.nurbs,...
                                            designParams,options.figBounds,2)
                                
        fname = [options.filePrefix,'_design_params_and_bounds'];
        set(fig, 'PaperPositionMode', 'auto');
        print('-djpeg ', '-r300', fname);
        saveas(fig, fname, 'fig');      
    else
        figure
        plot_channels_and_design_parameters(channels.contvty,channels.nurbs,...
                                            designParams,[],2)
    end
end
                              

varargout{1} = channels;
varargout{2} = bodySource;
varargout{3} = designParams;
varargout{4} = restrictedParams;

end
