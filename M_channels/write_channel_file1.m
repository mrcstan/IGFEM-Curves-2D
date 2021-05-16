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
% Omit this from channel input file{
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
% }
% diameters, <nch>
%
% channels.nurbs control points, 1, <nd>, <ncpts>
% xc_11 ... xc_nd1
% xc_12 ... xc_nd2
% ...
% xc_1ncpts ... xc_nd_ncpts
% channels.nurbs knot vector, 1, <nknots>
% xi ... xi_knots
% channels.contvty, 1
% i j 
% ...
% ...
% channels.nurbs control points, <nch>, <nd>, <ncpts>
% xc_11 ... xc_nd1
% xc_12 ... xc_nd2
% ...
% xc_1ncpts ... xc_nd_ncpts
% channels.nurbs knot vector, <nch>, <nknots>
% xi1 ... xi_knots
% channels.contvty, 1
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
function [] = write_channel_file(fname,channels,designParams,permission,headers)

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
fprintf(fileID,'end point number and temperature, %i\n', nInTemp);
for i = 1:nInTemp
    fprintf(fileID,'%i %g\n',channels.pt_temp(i,:));
end
if (channels.model == 1)
    modelType = 'mean temperature';
elseif (channels.model == 2)
    modelType = 'constant heat';
else
    error('unrecognized channel model')
end
fprintf(fileID,'model, %s\n',modelType);
fprintf(fileID,'mass in, %i, %g\n',channels.inletEndPoint,channels.massin);
fprintf(fileID,'heat capacity, %g\n',channels.heatCapacity);
fprintf(fileID,'viscosity, %g\n',channels.viscosity);
fprintf(fileID,'density, %g\n',channels.density);
fprintf(fileID,'pressure out, %i, %g\n',channels.pressureOutletEndPoint,...
                                        channels.pressureOut);
fprintf(fileID,'cross section, %s \n',channels.crossSection);

fprintf(fileID,'diameters, %i \n',nChannels);
fprintf(fileID,'%g\n',channels.diams);
for i = 1:nChannels
    if (any(channels.nurbs(i).coefs(4,:) ~= 1))
        dim = 4;
        formatStr = '%g %g %g %g\n';
    elseif (any(channels.nurbs(i).coefs(3,:)))
        dim = 3;
        formatStr = '%g %g %g\n';
    else
        dim = 2;
        formatStr = '%g %g \n';
    end
    fprintf(fileID,'nurbs control points, %i, %i, %i\n',i,dim,channels.nurbs(i).number);
    
    for j = 1:channels.nurbs(i).number
        fprintf(fileID,formatStr,channels.nurbs(i).coefs(1:dim,j));
    end
    
    fprintf(fileID,'nurbs knot vector, %i, %i\n', i, numel(channels.nurbs(i).knots));
    fprintf(fileID,'%g ',channels.nurbs(i).knots);
    fprintf(fileID,'\nconnectivity, %i\n', i);
    fprintf(fileID,'%i %i\n',channels.contvty(i,:));
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
        fprintf(fileID,'design, %i \n', i);
  
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