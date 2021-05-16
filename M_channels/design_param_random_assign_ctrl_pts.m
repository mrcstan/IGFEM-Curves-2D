%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Marcus Tan on 4/11/2015
%%% Last modified date: 4/11/2015
%%% Copyright 2015 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this function generates random initial designs by first generating as many
% random points within the large bounding box as the movable control points and
% then randomly assigns these random points as the control points of the
% channels
% ASSUMPTION:
%   designParams is ordered such that entries 2*i-1 and 2*i correspond to
%   the coordinates of the same control points
function sample = design_param_random_assign_ctrl_pts(inputFile,...
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
                                    

% randomize the initial values                

del = zeros(designParams.nParams+restrictedParams.nParams,1);
curSampleSize = 0;
sample = nan(designParams.nParams,nSamples);
if (isempty(nlconfun))
    warning('no nonlinear constraint function specified')
end
if (mod(designParams.nParams,2) ~= 0)
    error('function does not work if number of design parameters is not even')
else
    nParams2 = designParams.nParams/2;
end
tol = 1e-8;
while(curSampleSize < nSamples)
    fprintf('current sample size %i \n',curSampleSize)
    iniVals = randomized_bounded_values(designParams.bounds);
    for i = 1:maxTrials
        randInd = randperm(nParams2);
        randInd = [2*randInd-1;2*randInd];
        designParams.iniVals = iniVals(randInd(:));
        restrictedParams.iniVals ...
            = update_restricted_params_ini_vals(designParams, ...
                                                restrictedParams);
        channels = update_channels([designParams.iniVals;restrictedParams.iniVals],...
                                   designParams, ...
                                   restrictedParams, ...
                                   channels, ...
                                   'replace', ...
                                   false);
        if(channels_self_intersections(channels,tol))
            continue
        end
        if (isempty(nlconfun))
            curSampleSize = curSampleSize+1;
            sample(:,curSampleSize) = designParams.iniVals;
            break
        else
            fprintf('checking whether nonlinear constraints are satisfied at trial %i\n',i)
            try 
                if (any(nlconfun(del,designParams,restrictedParams,channels,nlcon) > 0))
                    continue
                else                                          
                    curSampleSize = curSampleSize+1;
                    sample(:,curSampleSize) = designParams.iniVals;
                    break
                end
            catch err
                 warning('generated samples may not satisfy nonlinear constraint function')
                 warning(getReport(err));
                 curSampleSize = curSampleSize+1;
                 sample(:,curSampleSize) = designParams.iniVals;
                 break
            end
        end
    end
    
end
   
end
