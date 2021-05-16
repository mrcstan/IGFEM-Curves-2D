%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Marcus Tan on 10/11/2014
%%% Last modified date: 2/15/2015
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
function sample = design_param_random_sampling(inputFile,...
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
                                    
generated = false;
% randomize the initial values                

del = zeros(designParams.nParams+restrictedParams.nParams,1);
curSampleSize = 0;
sample = nan(designParams.nParams,nSamples);
if (isempty(nlconfun))
    warning('no nonlinear constraint function specified')
end

for i = 1:maxTrials
    fprintf('trial %i \n',i)
     fprintf('current sample size %i \n',curSampleSize)
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
    if (isempty(nlconfun))
        curSampleSize = curSampleSize+1;
        sample(:,curSampleSize) = designParams.iniVals;
    else
        fprintf('checking whether nonlinear constraints are satisfied\n')
        try 
            if (any(nlconfun(del,designParams,restrictedParams,channels,nlcon,[],[],[]) > 0))
                continue
            else                                          
                curSampleSize = curSampleSize+1;
                sample(:,curSampleSize) = designParams.iniVals;
            end
        catch err
             warning('generated samples may not satisfy nonlinear constraint function')
             warning(getReport(err));
             curSampleSize = curSampleSize+1;
             sample(:,curSampleSize) = designParams.iniVals;
        end
    end
    if (curSampleSize == nSamples)
        generated = true;
        fprintf('number of attempts to get random channel = %i \n',i)
        break
    end
    
end
   


if (~generated)
    warning('maximum number of attempts to randomize channel has been exceeded')
end                       

end
