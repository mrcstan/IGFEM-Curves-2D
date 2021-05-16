%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Marcus Tan on 4/21/2015
%%% Copyright 2015 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this function creates a random parallel network with number of branches
% determined by looking at the inputFile
% Note: a keyword and a number is separated by a comma 
%       Other than the newline character or a blank, no other separator is
%       used
% INPUT: 
%  inputFile:
function sample = design_param_subgraph_monomorphism(inputFile,...
                                                     polygonFile,...
                                                     nSamples, ...
                                                     maxTrials, ...
                                                     coordBounds, ...
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

[vertex2channel,smallEdgeList,nOrgVertices,nVertices] ...
    = channel2vertex_vertex2channel_matrix(channels,'v2c');

% find original vertex coordinates
channel2vertex = channel2vertex_vertex2channel_matrix(channels,'c2v');
orgVertexCoord = nan(nVertices,2);
orgVertexCoord(1:nOrgVertices,:) = channels.pts;
for i= nOrgVertices+1:nVertices
    chanNum = find(channel2vertex(:,i));
    ctrlPtNum = channel2vertex(chanNum,i);
    orgVertexCoord(i,:) = channels.nurbs(chanNum).coefs(1:2,ctrlPtNum)';
end

if (isempty(nlconfun))
    warning('no nonlinear constraint function specified')
end
delx = coordBounds(1,2) - coordBounds(1,1);
dely = coordBounds(2,2) - coordBounds(2,1);
for attempts = 1:maxTrials
    fprintf('trial %i \n',attempts)
    fprintf('current sample size %i \n',curSampleSize)
    
    % ASSUME that the vertices 1 and nOrgVertices are the inlet and
    % outlet, respectively   
    vertexCoord = [channels.pts(1,:); ...
                   [delx*rand(nVertices-2,1)+coordBounds(1,1), ...
                    dely*rand(nVertices-2,1)+coordBounds(2,1)];...
                   channels.pts(nOrgVertices,:)];
    vertexCoord([nOrgVertices,nVertices],:) = vertexCoord([nVertices,nOrgVertices],:); 
    DT = delaunayTriangulation(vertexCoord);
    largeEdgeList = edges(DT);
    
    [small2largeCol1,small2largeCol2,matchNumber] ...
        = mx_subgraph_monomorphism(smallEdgeList',nVertices, ...
                                   largeEdgeList',nVertices, ...
                                   [1,nOrgVertices]);
     
    if (matchNumber <= 0)
        continue
    end
    % swap coordinates of the original channel network with the randomized
    % channel network
    vertexCoord(small2largeCol1,:) = vertexCoord(small2largeCol2,:); 
    % check if the orientation (ccw/cw) of the vertices are the same as before
    cycle = grCycleBasis(smallEdgeList);
    cycle1Edges = smallEdgeList(logical(cycle(:,1)),:); 
    G = sparse(cycle1Edges(:,1),cycle1Edges(:,2),1,nVertices,nVertices);
    discVertices = graphtraverse(G,cycle1Edges(1,1));
    if (ispolycw(orgVertexCoord(discVertices,1),orgVertexCoord(discVertices,2)) ...
            ~= ispolycw(vertexCoord(discVertices,1),vertexCoord(discVertices,2))) % different orientation (ccw/cw) of vertices 
        continue 
    end
    % update channel 
    for i = 1:nChannels
        nzInd = find(vertex2channel(:,i));
        channels.nurbs(i).coefs([1,2],vertex2channel(nzInd,i)) = vertexCoord(nzInd,:)';
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
