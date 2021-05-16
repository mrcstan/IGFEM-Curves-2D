%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Marcus Tan on 11/30/2014
%%% Last modified date: 2/14/2015
%%% Copyright 2014 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUT:
%   paramPairs: first column is the designParam
%               numbers and second column is the restricted param numbers
function [polygons,vertexCoords,vertices2params,...
          nRestrictedParams,rParamIniVals,paramPairs,sideTriangles] ...
                = channel_polygons(channels,nDesignParams,options)
if(nargin < 2)
    options.polygonFileName = [];
    options.mergeTriangles = true;
end

if (~isfield(options,'mergeTriangles'))
    options.mergeTriangles= false;
end

if (~isfield(options,'orientation'))
    options.orientation = 'NE';
end

nChannels = numel(channels.nurbs);
tot = sum(cat(1,channels.nurbs.number));
nInteriorPts = tot - 2*nChannels; 
vertexCoords = [channels.pts;
                nan(nInteriorPts,2)];
vertex_chanNum = nan(size(vertexCoords,1),1);
vertex_ctrlPtNum = nan(size(vertexCoords,1),1);
constraints = nan(2,tot - nChannels);
constraintCount = 0;
count = size(channels.pts,1);
for i = 1:nChannels
   vertex_chanNum(channels.contvty(i,1)) = i;
   vertex_chanNum(channels.contvty(i,2)) = i;
   vertex_ctrlPtNum(channels.contvty(i,1)) = 1;
   vertex_ctrlPtNum(channels.contvty(i,2)) = channels.nurbs(i).number;
   if (channels.nurbs(i).number > 2)
       newCount = count + channels.nurbs(i).number - 2;
       vertexCoords((count+1):newCount,:) = channels.nurbs(i).coefs(1:2,2:end-1)';
       vertex_chanNum((count+1):newCount) = i;
       vertex_ctrlPtNum((count+1):newCount) = 2:channels.nurbs(i).number-1;
       constraintCount = constraintCount + 1;
       constraints(:,constraintCount) = [channels.contvty(i,1);count+1];
       constraintCount = constraintCount + 1;
       constraints(:,constraintCount) = [newCount;channels.contvty(i,2)];
       for j = 1:(channels.nurbs(i).number-3)
           constraintCount = constraintCount + 1;
           constraints(:,constraintCount) = [count+j;count+j+1];
       end
       count = newCount;
   else
       constraintCount = constraintCount + 1;
       constraints(:,constraintCount) = channels.contvty(i,:)';
   end
end
constraints(:,constraintCount+1:end) = [];

% map vertices to design parameters
vertices2params = nan(size(vertexCoords));
for i = 1:size(vertexCoords,1)
    vertices2params(i,:) = channels.designParamNum{vertex_chanNum(i)}(1:2,vertex_ctrlPtNum(i))';     
end

% add additional vertices at the sides to provide additional support to the
% network and prevent intersections
ymin = min(vertexCoords(:,2));
ymax = max(vertexCoords(:,2));
nOldVertices = size(vertexCoords,1);
nBinEdges = 2*nOldVertices;
binEdges = linspace(ymin,ymax,nBinEdges);
[~,binInd] = histc(vertexCoords(:,2),binEdges);
[~,ia] = unique(binInd);
uniqueYs = vertexCoords(ia,2);
nUniqueYs = numel(uniqueYs);
xmin = min(vertexCoords(:,1));
xmax = max(vertexCoords(:,1));
vertexCoords = [vertexCoords;nan(2*(nUniqueYs-1),2)];

% create new left vertices
for i = 1:(nUniqueYs-1)
    vertexCoords(nOldVertices + i,:) = [xmin,uniqueYs(i)];
end

% create new right vertices
for i = 1:(nUniqueYs-1)
    j = nOldVertices + i + nUniqueYs - 1;
    vertexCoords(j,:) = [xmax,uniqueYs(i+1)];
end

%constraints = channels.contvty;
DT = delaunayTriangulation(vertexCoords,constraints');
TR = DT.ConnectivityList;

[elem2vertex,vertex2vertex] = adjacency_matrix(TR,size(vertexCoords,1));
paramPairs = nan(nDesignParams,2);
rParamIniVals = nan(nDesignParams,1);
nRestrictedParams = 0;
vertices2params = [vertices2params;nan(2*(nUniqueYs-1),2)];
tol = 1e-10;

leftVertices = find(abs(vertexCoords(:,1) - xmin) < tol);
for i = 1:numel(leftVertices)
    if (abs(vertexCoords(leftVertices(i),2) - ymax) < tol)
        continue
    end
    neighbourVertices = setdiff(find(vertex2vertex(:,leftVertices(i))),leftVertices(i));
    ind = abs(vertexCoords(neighbourVertices,2) - vertexCoords(leftVertices(i),2)) < tol;
    nRestrictedParams = nRestrictedParams + 1;
    nParams = nDesignParams + nRestrictedParams;
    vertices2params(leftVertices(i),2) = nParams;
    paramPairs(nRestrictedParams,:) = [vertices2params(neighbourVertices(ind),2),nParams];
    rParamIniVals(nRestrictedParams) = vertexCoords(leftVertices(i),2);
end

rightVertices = find(abs(vertexCoords(:,1) - xmax) < tol);
for i = 1:numel(rightVertices)
    if (abs(vertexCoords(rightVertices(i),2) - ymin) < tol)
        continue
    end
    neighbourVertices = setdiff(find(vertex2vertex(:,rightVertices(i))),rightVertices(i));
    ind = abs(vertexCoords(neighbourVertices,2) - vertexCoords(rightVertices(i),2)) < tol;
    nRestrictedParams = nRestrictedParams + 1;
    nParams = nDesignParams + nRestrictedParams;
    vertices2params(rightVertices(i),2) = nParams;  
    paramPairs(nRestrictedParams,:) = [vertices2params(neighbourVertices(ind),2),nParams];
    rParamIniVals(nRestrictedParams) = vertexCoords(rightVertices(i),2);
end
paramPairs(nRestrictedParams+1:end,:) = [];
rParamIniVals(nRestrictedParams+1:end) = [];

% ASSUMPTION: first point of channels.pts is the inlet and
%             last point of channels.pts is the outlet
%{
delTR = false(size(TR,1),1);
auxSet1 = [1,size(channels.pts,1)];
%auxSet2 = [1,size(channels.contvty,1)];
%auxSet3 = [2,1;1,2];
for j = 1:numel(auxSet1)
    indTR = vertexAttachments(DT,auxSet1(j));
    indTR = indTR{1};
    
    % remove all triangles associated with the inlet and outlet
    delTR(indTR) = true;
    
    % remove all but one triangle associated with the inlet and outlet
    %{
    if (numel(indTR) > 1)
        vec12 = channels.pts(channels.contvty(auxSet2(j),auxSet3(j,1)),:)...
               -channels.pts(channels.contvty(auxSet2(j),auxSet3(j,2)),:);
        vec12 = [vec12,0];
        leastNeg = -inf;
        leastNegInd = 1;
        for i = 1:numel(indTR)
            Xc = mean(vertexCoords(TR(indTR(i),:),:),1);
            vec1c = [Xc - channels.pts(channels.contvty(auxSet2(j),auxSet3(j,2)),:),0];
            crossProduct = cross(vec12,vec1c);
            component3 = crossProduct(3);
            if (component3 < 0 && component3 > leastNeg)
                leastNeg = component3;
                leastNegInd = i;
            end           
        end
        delTR(setdiff(indTR,indTR(leastNegInd))) = true;
    end
    %}
end


TR(delTR,:) = [];
DT = triangulation(TR,vertexCoords(:,1),vertexCoords(:,2));


tol = 1e-14;
leftVertices = vertexCoords(:,1) < mean(vertexCoords(:,1));
rightVertices = find(~leftVertices);
leftVertices = find(leftVertices);

[P,Q] = meshgrid(leftVertices,leftVertices);
mask = true(size(P));
mask = tril(mask,-1);
leftVertexPairs = [P(mask),Q(mask)];
mask = isConnected(DT,leftVertexPairs);
leftVertexPairs(mask,:) = [];

[P,Q] = meshgrid(rightVertices,rightVertices);
mask = true(size(P));
mask = tril(mask,-1);
rightVertexPairs = [P(mask),Q(mask)];
mask = isConnected(DT,rightVertexPairs);
rightVertexPairs(mask,:) = [];

topVertices = find(abs(vertexCoords(:,2) - max(vertexCoords(:,2))) < tol);
[P,Q] = meshgrid(topVertices,topVertices);
mask = true(size(P));
mask = tril(mask,-1);
topVertexPairs = [P(mask),Q(mask)];
mask = isConnected(DT,topVertexPairs);
topVertexPairs(mask,:) = [];

botVertices = find(abs(vertexCoords(:,2) - min(vertexCoords(:,2))) < tol);
[P,Q] = meshgrid(botVertices,botVertices);
mask = true(size(P));
mask = tril(mask,-1);
botVertexPairs = [P(mask),Q(mask)];
mask = isConnected(DT,botVertexPairs);
botVertexPairs(mask,:) = [];

vertexPairs = [leftVertexPairs;rightVertexPairs;
               topVertexPairs;botVertexPairs];
%}

if (options.mergeTriangles)
    [TR,rec] = merge_triangles(TR,vertexCoords,constraints);
    nRec = size(rec,1);
    nTR = size(TR,1);
    warning('sideTriangles is not valid for options.mergeTriangles = true')
    sideTriangles = [];
else
    nRec = 0;
    nTR = size(TR,1);
    
    % find side triangles 
    [rows,~] = find(elem2vertex(:,[leftVertices,rightVertices]));
    sideTriangles = unique(rows);
    %
    % remove triangles that contains
    % the inlet or outlet channels
    if (strcmpi(options.orientation,'NW'))
        delInd = false(numel(sideTriangles),1);
        for i = 1:numel(sideTriangles)
            if ((any(abs(vertexCoords(TR(sideTriangles(i),:),1) - xmin) < tol) ...
                 && nnz(abs(vertexCoords(TR(sideTriangles(i),:),2) - ymax) < tol) == 2) ...
                 || (any(abs(vertexCoords(TR(sideTriangles(i),:),1) - xmax) < tol) ...
                    &&  nnz(abs(vertexCoords(TR(sideTriangles(i),:),2) - ymin) < tol) == 2))
                delInd(i) = true;
            end
        end
        sideTriangles(delInd) = [];
    % swap orientation of triangles from NW to NE    
    elseif (strcmpi(options.orientation,'NE'))
        edgeNodes = find_edge_node(TR);
        dualedge = elems_sharing_edge_nodes(TR,size(vertexCoords,1)); 
        for i = 1:size(edgeNodes,1)
            delXX = vertexCoords(edgeNodes(i,2),:) - vertexCoords(edgeNodes(i,1),:);
            delXX = delXX/norm(delXX);
            if (abs(delXX(1)) > tol && abs(delXX(2)) > tol)
             sharedElems = [dualedge(edgeNodes(i,1),edgeNodes(i,2)),...
                            dualedge(edgeNodes(i,2),edgeNodes(i,1))];
             [edgeNodes(i,:),TR(sharedElems,:)]...
                 = flip_edge(edgeNodes(i,:),TR(sharedElems,:));
            end
        end
    else
        error('unknown orientation of triangles')
    end
end




polygons = struct('nVertices',0,'connectivity',cell(nTR+nRec,1));
for i = 1:nTR
    polygons(i).nVertices = 3;
    polygons(i).connectivity = TR(i,:);
end

for i = 1:nRec
    ii = i+nTR;
    polygons(ii).nVertices = 4;
    polygons(ii).connectivity = rec(i,:);
end


end