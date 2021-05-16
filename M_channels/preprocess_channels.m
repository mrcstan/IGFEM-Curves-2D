%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Marcus Tan on 10/11/2014
%%% Last modified date: 2/6/2015
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
%  options.sampleFile: Latin hypercube sampling file. if both boundsFile and
%                   sampleFile are present, sampleFile is used.
%  options.maxTrials: max number of trials to generate random configuration
%                     that satisfies nonlinear constraints
%  options.nlconfun: a function handle to the nonlinear constraint
%  options.nlcon: other inputs for nonlinear constraint including the rhs
%  options.polygonGenFile: file with same format as channel file to generate 
%                       polygons that are needed by some
%                       nonlinear constraints. this is guaranteed to work 
%                       only for reference configuration
%  options.polygonFig: polygons output file name
% OUTPUT:
%  channels:
%  designParams:
%  restrictedParams:
function varargout = preprocess_channels(inputFile,options)                                   
if (~ischar(inputFile))
    error('input file name must be a character array');
end

if (nargin < 2)
    options.figBounds = [];
    options.filePrefix = [];
    options.guimode = false;
    options.boundsFile = [];
    options.sampleFile = [];
    options.polygonFile = []; 
    options.maxTrials = 50;
    options.nlconfun = [];
    options.nlcon = [];
    
    options.polygonFig = [];
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
if (~isfield(options,'sampleFile'))
    options.sampleFile = [];
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
if (~isfield(options,'polygonFile'))
    options.polygonFile = [];
end

if (~isfield(options,'polygonFig'))
    options.polygonFig = [];
end

if (~isfield(options,'suppressDesignParams'))
    options.suppressDesignParams = false;
end

if (ischar(options.sampleFile) && ischar(options.boundsFile))
    warning('Latin hypercube sampling file is used instead of the bounds file')
    options.boundsFile = [];
end

if (ischar(options.sampleFile))
    if(~isfield(options,'sampleNum'))
        error('must specify sample number to select from Latin hypercube sampling file')
    end
    if(~isnumeric(options.sampleNum))
        error('sample number must be a numeric')
    end
end

[channels,designParams] = read_channels(inputFile);

if options.suppressDesignParams
    designParams.nParams = 0;
end
% design parameters
if (~options.guimode)
    [designParams,channels.designParamNum] ...
        = design_params2channel_params_map(designParams,...
                                           channels);
end


% randomize the initial values 
if (ischar(options.sampleFile))
    designParams.iniVals = read_sample_file(options.sampleFile,options.sampleNum);
    %restrictedParams.iniVals ...
    %    = update_restricted_params_ini_vals(designParams, ...
    %                                        restrictedParams);
end

if (ischar(options.polygonFile))
    [channels.polygons, channels.vertexCoords,...
     designParams.vertices2params,...
     restrictedParams.nParams, ...
     restrictedParams.iniVals, ...
     restrictedParams.paramPairs, ...
     sideTriangles] ...
        = read_polygon_file(options.polygonFile);
    
    [channels.polygons.isSideTriangle] = deal(false);
    [channels.polygons(sideTriangles).isSideTriangle] = deal(true);
    % make sure that the restricted params have the same values as the
    % corresponding params. 
    restrictedParams.iniVals ...
        = update_restricted_params_ini_vals(designParams, ...
                                            restrictedParams);
                                   
else
    restrictedParams.nParams = 0;
    restrictedParams.iniVals = [];
end

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
% need to calculate pressure, mass and mcf even for randomFlag == true
% because the update_channels function does not calculate those when those
% output are not requested
if (~options.guimode)
    pressure = [];
    mass = [];
    if (designParams.nParams)
        [channels,pressure,mass] ...
                = update_channels([designParams.iniVals;restrictedParams.iniVals], ...
                                  designParams, ...
                                  restrictedParams, ...
                                  channels, ...
                                  'replace', ...
                                  false);                   
        channels.Pin = pressure(channels.inletEndPoint);                      
         if (ischar(options.polygonFig))
             %fig = figure('visible','off');
             fig = figure;
             plot_channel_polygons(channels.polygons,...
                                   channels.vertexCoords,...
                                   channels.nurbs,...
                                   designParams,...
                                   false)
             saveas(fig, options.polygonFig, 'jpg');    
         end
     
    else
        ngpts = 10;
        [channels.length,channels.lineSegs,channels.segLengths]...
                    = nurbs_channel_lengths(channels,ngpts);    
        channels.vol = channel_volume(channels.length,...
                                      channels.diams,...
                                      channels.heights,...
                                      channels.crossSection);
        channels.branch_kinks = branching_points_n_kinks(channels);
        %channels.designParamNum = [];
        % only calculates the flow rate and pressure when channels.massin
        % is specified but not channels.mcf
        if ~isfield(channels,'mcf')
            [pressure,mass] = network_pressure_mass_flow_rate(channels.contvty,...
                                                              channels.nurbs,...
                                                              channels.diams,...
                                                              channels.heights,...
                                                              channels.viscosity,...
                                                              channels.inletEndPoint,...
                                                              channels.massin,...
                                                              channels.powerXdensity,...
                                                              channels.pressureOutletEndPoint,...
                                                              channels.pressureOut,...
                                                              channels.crossSection);
            channels.mcf = mass*channels.heatCapacity;  
            channels.Pin = pressure(channels.inletEndPoint);
        end
    end
    channels.Pin = pressure(channels.inletEndPoint);
    
    if (any(isnan(pressure)) || any(isnan(mass)))
        error('some nodal pressures or channel flow rates are NaN')
    end
end

if ~isfield(options,'selfIntersectTol')
    options.selfIntersectTol = 1e-10;
end
if(channels_self_intersections(channels,options.selfIntersectTol))
    error('channels self intersect')
end   

% plot pressure and mass flow rate in terms of kPa and ml/min, respectively
if (~options.guimode)
    if (ischar(options.filePrefix))
        %fig = figure('visible','off');
        %set(fig,'ResizeFcn','set(gcf,''visible'',''on'')'); 
        fig = figure;
        plot_channel_network(channels.contvty,channels.nurbs,pressure/1000.0,mass*1e6*60/1065)
        fname = [options.filePrefix,'_ini_mass_pressure'];
        %set(fig, 'PaperPositionMode', 'auto');
        %print('-djpeg ', '-r300', fname);
        saveas(fig, fname, 'jpg');
        saveas(fig, fname, 'fig');
    else
        figure
        plot_channel_network(channels.contvty,channels.nurbs,pressure/1000.0,mass*1e6*60/1065)
    end
end

if (~options.guimode && designParams.nParams)
    if (ischar(options.filePrefix))
        %fig = figure('visible','off');
        %set(fig,'ResizeFcn','set(gcf,''visible'',''on'')'); 
        fig = figure;
        plotOptions.figBounds = options.figBounds;
        plot_channels_and_design_parameters(channels.contvty, ...
                                            channels.pts, ...
                                            channels.nurbs,...
                                            designParams,plotOptions)
                                
        fname = [options.filePrefix,'_design_params_and_bounds'];
        %set(fig, 'PaperPositionMode', 'auto');
        %print('-djpeg ', '-r300', fname);
        saveas(fig, fname, 'jpg');
        saveas(fig, fname, 'fig');    
  
    else
        figure;     
        plotOptions.figBounds = options.figBounds;
        plot_channels_and_design_parameters(channels.contvty, ...
                                            channels.pts, ...
                                            channels.nurbs,...
                                            designParams,plotOptions)
    end
end
                              

varargout{1} = channels;
varargout{2} = designParams;
varargout{3} = restrictedParams;

end
