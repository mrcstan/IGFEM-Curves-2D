%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Marcus Tan on 10/11/2014
%%% Last modified date: 2/15/2015
%%% Copyright 2014 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read channel file input that contains bounds on design parameters and
% use latin hypercube to sample the control points within the bounds and 
% check that all of them satisfy a nonlinear constraint nlconfun
% Note: a keyword and a number is separated by a comma 
%       Other than the newline character or a blank, no other separator is
%       used
% INPUT: (if design parameters are not specified in input file only first
%         variables is relevant)
%  inputFile:

function lhs = design_param_latin_hypercube_sampling(inputFile,...
                                                     polygonFile,...
                                                     nSamples, ...
                                                     maxTrials,...
                                                     nlconfun,...
                                                     nlcon)

if (~ischar(inputFile))
    error('input file name must be a character array');
end



[channels,designParams] = read_channels(inputFile);

[designParams,channels.designParamNum] ...
    = design_params2channel_params_map(designParams,...
                                       channels);

%{
options.polygonFileName = [];
[channels.polygons, channels.vertexCoords,...
 designParams.vertices2params,...
 restrictedParams.nParams, ...
 restrictedParams.iniVals, ...
 restrictedParams.paramPairs] ...
    = channel_polygons(channels,designParams.nParams,options);
%}
 [channels.polygons, channels.vertexCoords,...
     designParams.vertices2params,...
     restrictedParams.nParams, ...
     restrictedParams.iniVals, ...
     restrictedParams.paramPairs, ...
     sideTriangles] ...
        = read_polygon_file(polygonFile);
  [channels.polygons.isSideTriangle] = deal(false);
  [channels.polygons(sideTriangles).isSideTriangle] = deal(true);      
% make sure that the restricted params have the same values as the
% corresponding params. 
restrictedParams.iniVals ...
    = update_restricted_params_ini_vals(designParams, ...
                                        restrictedParams);

                                   
del = zeros(designParams.nParams+restrictedParams.nParams,1); 
if (isempty(nlconfun))
    warning('no nonlinear constraint function specified')
end
for i = 1:maxTrials
    fprintf('trial %i \n',i)
    lhs = lhsdesign(nSamples,designParams.nParams)';
    lhs = bsxfun(@times,1-lhs,designParams.bounds(:,1)) ...
         +bsxfun(@times,lhs,designParams.bounds(:,2));
    generated = true;
    if (isa(nlconfun,'function_handle'))
        for j = 1:nSamples
            designParams.iniVals = lhs(:,j);
            restrictedParams.iniVals ...
                = update_restricted_params_ini_vals(designParams, ...
                                                    restrictedParams);
            channels = update_channels([designParams.iniVals;restrictedParams.iniVals],...
                                       designParams, ...
                                       restrictedParams, ...
                                       channels, ...
                                       'replace', ...
                                       false);
            if (~isempty(nlconfun))
                try     
                    if (any(nlconfun(del,designParams,restrictedParams,channels,nlcon) > 0))
                        generated = false;
                        break

                    end       
                catch err
                    warning('generated samples may not satisfy nonlinear constraint function')
                    warning(getReport(err));
                    break
                end
            end
        end
    else
        fprintf('nonlinear constraint function unavailable \n')
    end
    if (generated)
        fprintf('number of attempts to get random channel = %i \n',i)
        break
    end
end

if (~generated)
    warning('maximum number of attempts has been exceeded without generating a sample that satisfy the nonlinear constraints')
end
                     

end
