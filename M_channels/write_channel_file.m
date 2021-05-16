%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Marcus Tan on 11/5/2014
%%% Last modified date: 11/5/2014
%%% Copyright 2014 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% write channel input file from channels.nurbs variables
% INPUT:
%   fname: file name
%   channels.nurbs: a vector of channels.nurbs structure variable
% OUTPUT:
%   status: 0: success
%           -1: failure
% Note: a keyword and a number is separated by a comma 
%       Other than the newline character or a blank, no other separator is
%       used
% FORMAT:
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
% xc_11 ... xc_nd1
% xc_12 ... xc_nd2
% ...
% xc_1ncpts ... xc_nd_ncpts
% nurbs knot vector, 1, <nknots>
% xi1 ... xi_knots
% connectivity, 1
% i j 
% ...
% ...
% nurbs control points
% xc_11 ... xc_nd1
% xc_12 ... xc_nd2
% ...
% xc_1ncpts ... xc_nd_ncpts
% nurbs knot vector
% xi1 ... xi_knots
% connectivity
% i j 
% ...
% ...
% <empty line>
% nurbs control points
% xc_11 ... xc_nd1
% xc_12 ... xc_nd2
% ...
% xc_1ncpts ... xc_nd_ncpts
% nurbs knot vector
% xi1 ... xi_knots
% connectivity
% i j 
%
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
function write_channel_file(fname,channels,designParams,permission,headers)

if (~ischar(fname))
    error('output file name must be a character array');
end
if (~strcmp(permission,'w') && ~strcmp(permission,'a'))
    warning('cannot write channel file')
    return
end

nChannels = numel(channels.nurbs);

if (size(channels.contvty,1) ~= nChannels || numel(channels.diams) ~= nChannels)
    error('size(channels.contvty,1) and numel(channels.diams) not equal numel(channels.nurbs)')
end

fileID = fopen(fname,permission);
if (fileID == -1)
    warning('cannot create or open file')
end
if (~isempty(headers))
    if (~iscell(headers))
        error('headers must a cell array of strings')
    end
    for i = 1:numel(headers)
        fprintf(fileID,'%s \n',headers{i});
    end 
end
fprintf(fileID,'number of channels, %i\n',nChannels);
nInTemp = size(channels.pt_temp,1);
fprintf(fileID,'\n');
fprintf(fileID,'end point number and temperature, %i\n', nInTemp);
for i = 1:nInTemp
    fprintf(fileID,'%i %g\n',channels.pt_temp(i,:));
end
fprintf(fileID,'\n');
if (channels.model == 1)
    modelType = 'mean temperature';
elseif (channels.model == 2)
    modelType = 'constant heat';
else
    error('unrecognized channel model')
end
fprintf(fileID,'model, %s\n',modelType);
if ~isempty(channels.massin)
    
    fprintf(fileID,'mass in, %i, %g\n',channels.inletEndPoint,channels.massin);
    if ~isfield(channels,'powerXdensity') ...
            || numel(channels.powerXdensity) ~= numel(channels.inletEndPoint)
        powerXdensity = [];
        powerXdensityInlet = [];
    else        
        powerXdensity = channels.powerXdensity;
        powerXdensityInlet = channels.inletEndPoint;
    end
    fprintf(fileID,'# powerXdensity, %i, %g\n',powerXdensityInlet,powerXdensity);
elseif isfield(channels,'powerXdensity')
    fprintf(fileID,'powerXdensity, %i, %g\n',channels.inletEndPoint,channels.powerXdensity);
    if isfield(channels,'massin') ...
            || numel(channels.massin) ~= numel(channels.inletEndPoint)
        massin = [];
        massInlet = [];
    else        
        massin = channels.massin;
        massInlet = channels.inletEndPoint;
    end
    fprintf(fileID,'# mass in, %i, %g\n',massInlet,massin);
end
fprintf(fileID,'heat capacity, %g\n',channels.heatCapacity);
fprintf(fileID,'viscosity, %g\n',channels.viscosity);
if isfield(channels,'viscosityModel')
    viscosityModel = channels.viscosityModel;
else
    viscosityModel = 1;
end
fprintf(fileID,'viscosity model, %g\n',viscosityModel);
fprintf(fileID,'density, %g\n',channels.density);
fprintf(fileID,'pressure out, %i, %g\n',channels.pressureOutletEndPoint,...
                                        channels.pressureOut);
fprintf(fileID,'cross section, %s \n',channels.crossSection);

fprintf(fileID,'\ndiameters \n');
fprintf(fileID,'%g\n',channels.diams);
fprintf(fileID,'\n');

if isfield(channels,'heights') && ~isempty(channels.heights)
    fprintf(fileID,'heights \n');
    fprintf(fileID,'%g\n',channels.heights);
    fprintf(fileID,'\n');
end

for i = 1:nChannels
    fprintf(fileID,'nurbs control points\n');
    
    for j = 1:channels.nurbs(i).number
        if (any(channels.nurbs(i).coefs(4,:) ~= 1))
            fprintf(fileID,'%g %g %g %g\n',[channels.nurbs(i).coefs(1:3,j)...
                                 ./channels.nurbs(i).coefs(4,j);...
                                  channels.nurbs(i).coefs(4,j)]);
        elseif (any(channels.nurbs(i).coefs(3,:)))
            fprintf(fileID,'%g %g %g\n',channels.nurbs(i).coefs(1:3,j));
        else
            fprintf(fileID,'%g %g\n',channels.nurbs(i).coefs(1:2,j));
        end
    end  
    fprintf(fileID,'nurbs knot vector\n');
    fprintf(fileID,'%g ',channels.nurbs(i).knots);
    fprintf(fileID,'\nconnectivity\n');
    fprintf(fileID,'%i %i\n',channels.contvty(i,:));
    fprintf(fileID,'\n');
end

if (~isempty(designParams))
    fprintf(fileID,'number of design variables, %i \n', designParams.nParams); 
    for i = 1:designParams.nParams
        if (isempty(designParams.channelNum{i}) ...
            && isempty(designParams.ctrlPtNum{i}) ...
            && isnan(designParams.ctrlPtnum(i)) ...
            && any(isinf(designParams.bounds(i,:))))
            continue
        end
        fprintf(fileID,'design\n');
  
        if (strcmpi(designParams.type{i},'CTRL_PT'))
            if (designParams.ctrlPtDim(i) == 1)
                str = 'channel, %i, cpt, %i, x, %g, %g \n';
            elseif (designParams.ctrlPtDim(i) == 2)
                str = 'channel, %i, cpt, %i, y, %g, %g \n';
            end
           
            fprintf(fileID,str,designParams.channelNum{i}(1),...
                        designParams.ctrlPtNum{i}(1),...
                        designParams.bounds(i,1),...
                        designParams.bounds(i,2));
        elseif (strcmpi(designParams.type{i},'DIAM'))
            str = 'channel, %i, diam, %g, %g \n';
            fprintf(fileID,str,designParams.channelNum{i}(1),...
                        designParams.bounds(i,1),...
                        designParams.bounds(i,2));
        else
            error('unknown design parameter type')
        end
        
    end
   
end
fclose(fileID);
end