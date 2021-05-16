%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Marcus Tan on 3/18/2015
%%% Copyright 2015 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this function creates a random parallel network with number of branches
% determined by looking at the inputFile
% Note: a keyword and a number is separated by a comma 
%       Other than the newline character or a blank, no other separator is
%       used
% INPUT: (if design parameters are not specified in input file only first
%         variables is relevant)
%  inputFile:
function sample = design_param_random_parallel_network(inputFile,...
                                                       polygonFile,...
                                                       nSamples, ...
                                                       maxTrials, ...
                                                       coordBounds, ...
                                                       lengthBounds, ...
                                                       minAngle,...
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
nChannels = size(channels.contvty,1);
nPts = max(channels.contvty(:));
if (mod(nPts,2) ~= 0)
    error('number of channel end points should be an even number')
else
    nBranches = nPts/2 - 1;
end
maxTrials1 = 100;

if (isempty(nlconfun))
    warning('no nonlinear constraint function specified')
end
for attempts = 1:maxTrials
    fprintf('trial %i \n',attempts)
    fprintf('current sample size %i \n',curSampleSize)
    
    [sortedContvty,sortedInd] = sortrows(channels.contvty);
    % ASSUME inlet is sortedContvty(1,1)
    % Also ASSUME outlet is sortedContvty(end,end)
    [vertices, ~ ,generated] ...
        = generate_random_parallel_network(channels.pts(1,:), ...
                                           channels.pts(end,:), ...
                                           nBranches, ...
                                           coordBounds, ...
                                           lengthBounds, ...
                                           minAngle, ...
                                           maxTrials1, ...
                                           false);
    if (generated)
        % first line segment
        %channels.pts(sortedContvty(1,2),:) = vertices(2,:);
        chan = sortedInd(1);
        channels.nurbs(chan).coefs(1:2,2) = vertices(2,:)';      
        
        % first loop
        for i = 2:3
            chan = sortedInd(i);
            if (channels.nurbs(chan).number > 2)
                channels.nurbs(chan).coefs(1:2,1) = vertices(2,:)';
                channels.nurbs(chan).coefs(1:2,2) = vertices(3,:)';
                channels.nurbs(chan).coefs(1:2,3) = vertices(5,:)';
                %channels.pts(sortedContvty(i,1),:) = vertices(2,:);
                %channels.pts(sortedContvty(i,2),:) = vertices(5,:);
            else
                channels.nurbs(chan).coefs(1:2,1) = vertices(2,:)';
                channels.nurbs(chan).coefs(1:2,2) = vertices(4,:)';
                %channels.pts(sortedContvty(i,1),:) = vertices(2,:);
                %channels.pts(sortedContvty(i,2),:) = vertices(4,:);
            end       
        end
        
        % middle loops
        for i = 4:(numel(sortedInd)- 3)
            chan = sortedInd(i);
            channels.nurbs(chan).coefs(1:2,1) = vertices(sortedContvty(i,1)+1,:)';
            channels.nurbs(chan).coefs(1:2,2) = vertices(sortedContvty(i,2)+1,:)';
            %channels.pts(sortedContvty(i,1),:) = vertices(sortedContvty(i,1)+1,:)';
            %channels.pts(sortedContvty(i,2),:) = vertices(sortedContvty(i,2)+1,:)';
        end
        
        % last loop
        for i = (numel(sortedInd)-2):(numel(sortedInd)-1)
            chan = sortedInd(i);
            if (channels.nurbs(chan).number > 2)
                channels.nurbs(chan).coefs(1:2,1) = vertices(sortedContvty(i,1)+1,:)';
                channels.nurbs(chan).coefs(1:2,2) = vertices(sortedContvty(i,1)+3,:)';
                channels.nurbs(chan).coefs(1:2,3) = vertices(sortedContvty(i,2)+2,:)';
                %channels.pts(sortedContvty(i,1),:) = vertices(sortedContvty(i,1)+1,:)';
                %channels.pts(sortedContvty(i,2),:) = vertices(sortedContvty(i,2)+4,:)';
            else
                channels.nurbs(chan).coefs(1:2,1) = vertices(sortedContvty(i,1)+1,:)';
                channels.nurbs(chan).coefs(1:2,2) = vertices(sortedContvty(i,2)+2,:)';
                %channels.pts(sortedContvty(i,1),:) = vertices(sortedContvty(i,1)+1,:)';
                %channels.pts(sortedContvty(i,2),:) = vertices(sortedContvty(i,2)+2,:)';
            end
        end
        i = numel(sortedInd);
        chan = sortedInd(i);
        channels.nurbs(chan).coefs(1:2,1) = vertices(sortedContvty(i,1)+2,:)';
    else
        continue
    end
    for i = 1:nChannels
        ind = ~isnan(channels.designParamNum{i});
        coefs = channels.nurbs(i).coefs([1,2,4],:);
        designParams.iniVals(channels.designParamNum{i}(ind)) = coefs(ind);
    end
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
            if (any(nlconfun(del,designParams,restrictedParams,channels,nlcon) > 0))
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
        fprintf('number of attempts to get random channel = %i \n',attempts)
        break
    end
    
end
   


if (~generated)
    warning('maximum number of attempts to randomize channel has been exceeded')
end                       

end
