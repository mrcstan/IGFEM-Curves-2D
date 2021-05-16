%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Marcus Tan on 8/25/2013
%%% Last modified date: 12/24/2014
%%% Copyright 2013 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% this function find elements with enriched nodes.
% ASSUMPTIONS: i) parent element is triangular
% LIMITATIONS: i) cannot handle triple intersection of a channel with an
%                 element, with one intersection per edge and one
%                 intersection point not connected with the rest
% INPUT:
% parent(i).edge: a length parent.nNodes array of edges upon which the nodes
%                 lie. For nodes at the vertices, the convention of the edge
%                 between vertex 1 and 2 is assigned to vertex 1, the edge
%                 between vertex 2 and 3 being assigned to vertex 2 and the
%                 edge between vertex 3 and 1 being assigned to vertex 3 is
%                 applied.
% parent(i).locEdge:
% parent(i).seg: a length parent.nNodes array of cells of the nurbs segments
%                passing through the nodes of the corresponding indices. seg
%                is an empty cell if a node is not an intersection point
% parent(i).nurbsParam: a length parent.nNodes array of cells of the nurbs parametric 
%                    coordinate of the intersection point.
% parent(i).coincident: boolean vector indicating if an edge has a
%                       coindicent channel
% elem_nodes
% elem_edges
% nodeCoords
% nOriginalNodes
% junc.nJunc
% junc.node
% junc.seg
% channels.nurbs
% INPUT-OUTPUT:
% parent(i).nSeg: a length parent.nNodes array of the number of nurbs curve 
%                 segments passing through the nodes of the corresponding
%                 indices. nSeg is 0 if a node is not an intersection point.
%                 output for the purpose of plottting in
%                 plot_mesh_nurbs_parent
% OUTPUT:
% 1 <= i <= number of elements
% parent(i).type:  0 = regular element without channel
%                  1 = regular element with channel
%                  2 = polynomial IGFEM element with channel
%                  3 = NURBS IGFEM element with channel
% parent(i).hasJunc: true if it contains a branching point
%                    When parent.isSpecial == true or parent.hasJunc == true, 
%                    the triangular child elements are not
%                    merged to form quadrilateral child elements whenever
%                    possible. (REMOVED)
% parent(i).isSpecial: an element is a special element if
%                   a single line source coincide with 2 edges of an
%                   element (REMOVED)
% parent(i).multItrsect: true if there is more than 2 intersections per edge
%                     and at least one intersection is not an original
%                     node. if an edge belongs to an element with branching, 
%                     this flag is false even if the edge has multiple intersections.
% parent(i).nItrsect: number of intersection points
% parent(i).nodes: 1D array of the global node numbers (both enrichment nodes 
%                  and original nodes) of the parent element. 
%                  The nodes on the edges are arranged in counter-clockwise order
%                  Interior nodes are listed last in no particular order.

% parent(i).nurbsSeg: a length parent.nNodes array of cells of the nurbs data 
%                  structure of the curve segments
% parent(i).dnurbsSeg:
% parent(i).channelNum: a row vector of channel numbers
% parent(i).channelNodes: 2 x nChannels array of the channels in the parent
%                         element where the elements of the array are terminal 
%                         nodes of the channels. The nodes specified using the 
%                         global numbers.
%                         The parametric coordinate of the nurbs
%                         segment increases from the first column to the second
%                         column. 
% parent(i).channelLocNodes: 2 x nChannels array of the channels in the parent
%                            element where the elements of the array are terminal 
%                            nodes of the channels. The nodes specified using the
%                            local number with respect to the parent element. 
%                            The parametric coordinate of the nurbs
%                            segment increases from the first column to the second
%                            column. 
% parent(i).channelNurbsParam: 2 x nChannels matrix of the channels in the parent
%                               element where the elements of the array are terminal 
%                               nurbs parameters of the channels. 
%                               The parametric coordinate of the nurbs segment 
%                               increases from the first row to the second
%                               row. 
% parent(i).channelLocEdge:
% parent(i).collinearLoc: local nodes of numbers that are collinear
% parent(i).cstrLocNodes: local number of constraint nodes
% parent(i).cstrVals: constraint values
% parent(i).cstrRows: global equation numbers of additional equations for 
%                    Lagrange multiplier method. At this point in time,
%                    the total number of nodes is not final. Only after the
%                    child elements are constructed is the total number of
%                    nodes final and this final number is added to cstrRows
%                    to get the correct global equation numbers
function [parent,cstrElems,nIGFEMelems] ...
                = parent_elements_nurbs(parent, ...
                                        elem_nodes, ....
                                        elem_material, ...
                                        material,...
                                        elem_branch_kinks, ...
                                        nodeCoords, ...
                                        nOriginalNodes,...
                                        constraint, ...
                                        channels, ...
                                        nurbsParamTol, ...
                                        polyIGFEM,...
                                        calcVel)
if (calcVel)
     fields2remove = {'coincident','seg','locEdge','nurbsParam',...
                     'multItrsect','vx','vy','auxDesignParamNum'};
else
    fields2remove = {'coincident','seg','locEdge','nurbsParam',...
                     'multItrsect'};
end

% add fields to parent
[parent.nurbsSeg,...
 parent.dnurbsSeg, ...
 parent.channelNum, ...
 parent.channelNodes, ...
 parent.channelLocNodes, ...
 parent.channelNurbsParam, ...
 parent.channelLocEdge,...
 parent.collinearLoc,...
 parent.cstrLocNodes,...
 parent.cstrVals,...
 parent.cstrRows,...
 parent.conductivity] = deal([]);

cstrElems = zeros(constraint.number,1);
countCstrElems = 0;
nConstraints = numel(constraint.temp_node);

nIGFEMelems = 0;

% determine if candidate parent element is a parent element
for i = 1:size(elem_nodes,1) 
     parent(i).conductivity = material(elem_material(i)).conductivity; 
    if parent(i).nItrsect > 1
        %hasJunc = false; % indicates whether an element has a branching point or a kink
        % decide if element becomes a regular element with line source or a
        % parent element with child elements
        % candidate is a parent element because more than one of its edges
        % are cut, splitting the original element.
        if(any(parent(i).nodes > nOriginalNodes) ...
                || ~isempty(elem_branch_kinks(i).bkNums))
            if (polyIGFEM)
                parent(i).type = 2;
            else
                parent(i).type = 3;
            end
            nIGFEMelems = nIGFEMelems + 1;
        else
           
            parent(i).type = 1;
        end

        % add original nodes of element which are not already in the 
        % list  to parent element
        [parent(i).nodes,~,ib]...
            =union(parent(i).nodes,elem_nodes(i,:),'stable');
        nPad=length(ib);
        parent(i).nSeg=[parent(i).nSeg,zeros(1,nPad)];            
        parent(i).seg=[parent(i).seg,cell(1,nPad)];
        parent(i).nurbsParam=[parent(i).nurbsParam,cell(1,nPad)];
        %parent(i).edge=[parent(i).edge,elem_edges(i,ib)];
        parent(i).locEdge=[parent(i).locEdge,ib'];
        
        if (calcVel)
            parent(i).vx = [parent(i).vx,cell(1,nPad)];
            parent(i).vy = [parent(i).vy,cell(1,nPad)];
            parent(i).auxDesignParamNum = [parent(i).auxDesignParamNum,cell(1,nPad)];
        end
        % arrange the  nodes are in CCW order
        % exclude branching point/junc when arranging
        hullInd = my_convHull(nodeCoords(parent(i).nodes,1),...
                           nodeCoords(parent(i).nodes,2));

        parent(i).nodes = parent(i).nodes(hullInd);
        parent(i).nSeg = parent(i).nSeg(hullInd);
        parent(i).seg = parent(i).seg(hullInd);
        parent(i).nurbsParam = parent(i).nurbsParam(hullInd);
        % parent(i).edge = parent(i).edge(hullInd);
        parent(i).locEdge = parent(i).locEdge(hullInd);
        if (calcVel)
            parent(i).vx = parent(i).vx(hullInd);
            parent(i).vy = parent(i).vy(hullInd);
            parent(i).auxDesignParamNum = parent(i).auxDesignParamNum(hullInd);
        end
         % arrange the nodes such that the original nodes appear first
        [~,~,iPON] = intersect(elem_nodes(i,:),...
                               parent(i).nodes,'stable');
        [enrichNode,iEN] = setdiff(parent(i).nodes,elem_nodes(i,:),'stable');
        parent(i).nodes = [elem_nodes(i,:),enrichNode];
        parent(i).nSeg = [parent(i).nSeg(iPON), ...
                               parent(i).nSeg(iEN)];
        parent(i).seg = [parent(i).seg(iPON), ...
                              parent(i).seg(iEN)];
        parent(i).nurbsParam = [parent(i).nurbsParam(iPON),...
                                     parent(i).nurbsParam(iEN)];
        %parent(i).edge = [elem_edges(i,:), ...
        %                       parent(i).edge(iEN)];        
        parent(i).locEdge = [[1,2,3],parent(i).locEdge(iEN)];
        
        if (calcVel)
            parent(i).vx = [parent(i).vx(iPON),parent(i).vx(iEN)];
            parent(i).vy = [parent(i).vy(iPON),parent(i).vy(iEN)];
            parent(i).auxDesignParamNum = [parent(i).auxDesignParamNum(iPON),...
                                        parent(i).auxDesignParamNum(iEN)];
        end
       
        
        if ~isempty(elem_branch_kinks(i).bkNums)
            nBKs = numel(elem_branch_kinks(i).bkNums);
            kOld = numel(parent(i).nodes);
            parent(i).nItrsect = parent(i).nItrsect + nBKs;
            parent(i).nodes = [parent(i).nodes,elem_branch_kinks(i).nodes];
            parent(i).nSeg = [parent(i).nSeg,ones(1,nBKs)];
            parent(i).seg = [parent(i).seg,cell(1,nBKs)];
            parent(i).nurbsParam = [parent(i).nurbsParam, cell(1,nBKs)];
            % parent(i).edge = [parent(i).edge,zeros(1,nBKs)];
            parent(i).locEdge = [parent(i).locEdge,zeros(1,nBKs)];
            if calcVel
                parent(i).vx = [parent(i).vx,cell(1,2*nBKs)];
                parent(i).vy = [parent(i).vy,cell(1,2*nBKs)];
            end
            
            for k = 1:nBKs
                bkNum = elem_branch_kinks(i).bkNums(k);
                m = kOld + k;
                % kink
                chanNums = channels.branch_kinks.channelNums{bkNum};
                parent(i).nSeg(m) = numel(chanNums);
                parent(i).seg{m} = chanNums;
                parent(i).nurbsParam{m} = channels.branch_kinks.nurbsParams{bkNum}; 
                if calcVel    
                    parent(i).vx{m} = [1,0];
                    parent(i).vy{m} = [0,1];
                    parent(i).auxDesignParamNum{m} ...
                        = channels.designParamNum{chanNums(1)}...
                                (1:2,channels.branch_kinks.ctrlPtNums{bkNum}(1))';
                end                      
            end                                      
        end
        
        [parent(i).channelLocNodes,...
         parent(i).channelNodes,...
         parent(i).channelNurbsParam,...
         parent(i).channelNum]...
                = find_channels(parent(i),channels.nurbs,nodeCoords(elem_nodes(i,:),:));
        
        
        parent(i).channelLocEdge = channel_local_edge(parent(i));

        parent(i) = remove_non_coincident_channels(parent(i));
    
        %parent(i).collinearLoc...
        %        = find_collinear_sets(parent(i),elem_edges(i,:)); 



        % add constraint nodes 
        if (countCstrElems < nConstraints && parent(i).type > 1)
            [Lia,locNode] = ismember(constraint.temp_node,parent(i).nodes);
            if (any(Lia))
                parent(i).cstrLocNodes = locNode(locNode>0);
                parent(i).cstrVals = constraint.temp_value(Lia);
                %parent(i).nConstraint = length(parent(i).cstrLocNodes);
                % must add the total number of nodes after child elements
                % are constructed
                parent(i).cstrRows = find(Lia);
                countCstrElems = countCstrElems + 1;
                cstrElems(countCstrElems) = i;
            end
        end


        % assign nurbs curve segment to parent element
        if (~polyIGFEM)
            parent(i).nurbsSeg = channels.nurbs;
            parent(i).dnurbsSeg = channels.nurbs;
            nChannels = numel(parent(i).channelNum); 

            for j=1:nChannels            
                xi1 = parent(i).channelNurbsParam(1,j);
                xi2 = parent(i).channelNurbsParam(2,j);                
                % if there is more than 2 line sources, just use curve
                % fitting to find a 2nd order linear nurbs curve with 2 control
                % points.
                if((parent(i).multItrsect) ...
                        && channels.nurbs(parent(i).channelNum(j)).order > 2)
                    fprintf('replace curves with straight segments in parent element %i\n',i)
                    parent(i).nurbsSeg(j) ...
                            = nurbs_curve_fit_segment(...
                                    channels.nurbs(parent(i).channelNum(j)),...
                                    [xi1,xi2],2,2);
                else    

                    parent(i).nurbsSeg(j) ...
                                = subint_extractor(...
                                        channels.nurbs(parent(i).channelNum(j)),...
                                        xi1, xi2,nurbsParamTol);

                    %nurbsSeg=nurbs_curve_fit_segment(...
                    %    channels.nurbs(parent(i).channelNum(j)),[xi1,xi2],3,3);

                end

                parent(i).dnurbsSeg(j) = nrbderiv(parent(i).nurbsSeg(j));
            end          
            parent(i).nurbsSeg = parent(i).nurbsSeg(1:nChannels);
            parent(i).dnurbsSeg = parent(i).dnurbsSeg(1:nChannels);
        end
        
        if(calcVel)      
            parent(i).designParamNum = unique(cat(2,parent(i).auxDesignParamNum{:}));
            parent(i).designParamNum(isnan(parent(i).designParamNum)) = [];
            nUniqueParams = numel(parent(i).designParamNum);
            nEnrichNodes = numel(parent(i).nodes) - 3;
            parent(i).vel = zeros(2,nEnrichNodes,nUniqueParams);
            if (~isempty(parent(i).designParamNum))         
                for j = 4:numel(parent(i).nodes)
                    [isDefinedInd,paramLoc] = ismember(parent(i).auxDesignParamNum{j}, ...
                                                       parent(i).designParamNum);
                    paramLoc = paramLoc(paramLoc > 0);
                    parent(i).vel(1:2,j-3,paramLoc) = [parent(i).vx{j}(isDefinedInd);
                                                       parent(i).vy{j}(isDefinedInd)];
                end
                % remove design parameters that have all velocities == 0
                delParamInd = all(squeeze(all(parent(i).vel == 0,1)));
                if (~isempty(delParamInd))
                    parent(i).designParamNum(delParamInd) = [];
                    parent(i).vel(:,:,delParamInd) = [];
                end
            end      
        end
    end     
end
% remove fields that are no longer needed
parent = rmfield(parent,fields2remove);
end

% this function finds the end nodoes of line sources in a parent element
% and the label of the line sources
% OUTPUT: 
% channelLocNodes: a nx2 array of the local number of the end nodes of line 
%                sources. n is the number of line sources.
% channelNodes: a nx2 array of the global number of the end nodes of line
%               source.
% channelNurbsParam: nurbs parameter of the end points of the line source
% channelNum: a length n array of the line source segment numbers

function [channelLocNodes,channelNodes,channelNurbsParam,channelNum,isSpecial]...
          =find_channels(parent,channelNurbs,Xelem)
[uniqueSeg,~,ib] = unique(cat(2,parent.seg{:}));

nChannels = numel(uniqueSeg);
channelLocNodes = inf(2,nChannels);
channelNodes = inf(2,nChannels);
channelNurbsParam = inf(2,nChannels);
channelNum=uniqueSeg';
tempN=zeros(nChannels,1);
cumSumNSeg=cumsum(parent.nSeg);
for i = 1:numel(parent.nodes)
    if(parent.nSeg(i)>0)
        for j=1:parent.nSeg(i)
            if(i==1)
                k=j;
            else
                k=cumSumNSeg(i-1)+j;
            end
            tempN(ib(k)) = tempN(ib(k))+1;
            if (tempN(ib(k)) > size(channelNodes,1))
                channelLocNodes(tempN(ib(k)),:) = inf;
                channelNodes(tempN(ib(k)),:) = inf;
                channelNurbsParam(tempN(ib(k)),:) = inf;
            end
            channelLocNodes(tempN(ib(k)),ib(k)) = i;
            channelNodes(tempN(ib(k)),ib(k)) = parent.nodes(i);
            %if (isfield(parent,'nurbsParam'))
            channelNurbsParam(tempN(ib(k)),ib(k)) = parent.nurbsParam{i}(j);
            %end
        end
    end
end

nChannels=size(channelNodes,2);
% sort the channel nodes in order of ascending NURBS parameters
[channelNurbsParam,sortInd] = sort(channelNurbsParam);
for i = 1:nChannels
    channelLocNodes(:,i) = channelLocNodes(sortInd(:,i),i);
    channelNodes(:,i) = channelNodes(sortInd(:,i),i);
end
%{
sortInd = bsxfun(@plus,sortInd,(0:2:2*(nChannels-1)));
channelLocNodes = channelLocNodes(sortInd);
channelNodes = channelNodes(sortInd);
%}
% consider the case where a line source has three or more intersection points
nIntersect2=size(channelNodes,1)-2;
isSpecial = false; 
%
if(nIntersect2>0)    
    isSpecial = true;
    for i=1:nChannels
        %{
        % arrange intersection points in increasing nurbs parameter order
        xiParam = cell2mat(parent.nurbsParam(channelLocNodes(:,i)));
        [~,sortInd] = sort(xiParam);
        channelLocNodes(:,i) = channelLocNodes(sortInd,i); 
        channelNodes(:,i) = channelNodes(sortInd,i);
        channelNurbsParam(:,i) = channelNurbsParam(sortInd,i);
        %}
        for j=1:nIntersect2
            nChannels=nChannels+1;
            channelLocNodes(1:2,nChannels)=[channelLocNodes(1+j,i);channelLocNodes(2+j,i)];
            channelNodes(1:2,nChannels)=[channelNodes(1+j,i);channelNodes(2+j,i)];        
            channelNurbsParam(1:2,nChannels)=[channelNurbsParam(1+j,i);channelNurbsParam(2+j,i)];
            channelNum(nChannels)=channelNum(i);
        end
    end
end

% delete rows beyond 2nd row.
%
if(size(channelLocNodes,1)>2)
    channelLocNodes = channelLocNodes(1:2,:);
    channelNodes = channelNodes(1:2,:);
    channelNurbsParam = channelNurbsParam(1:2,:);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% delete cols with one zero or more.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
delcol = any(channelLocNodes == inf);
channelLocNodes(:,delcol) = [];
channelNodes(:,delcol) = [];
channelNum(delcol) = [];
channelNurbsParam(:,delcol) = [];

% check if channel segment between parent(i).channelNodes(1:2,j) lies within element
% if not delete parent(i).channelNodes(1:2,j)
% this is to prevent errornously ascribing a curved channel that has
% multiple nodes on the edges of an element but does not lie within the
% element
if (numel(unique(channelNum)) ~= numel(channelNum))
    nElemChannels = numel(channelNum); 
    delInd = false(nElemChannels,1);
    % here I took the shortcut and check all channels in the element 
    % instead of just the repeated channels 
    % I sample a point evaluated at  0.5*(channelNurbsParam(1,j)+channelNurbsParam(2,j))
    % and check if the point is in the element. if it is not, I consider
    % the channel segment to lie outside the element
    for j = 1:nElemChannels
        pt = nrbeval(channelNurbs(channelNum(j)), ...
                     0.5*(channelNurbsParam(1,j)+channelNurbsParam(2,j)));
        if (~inpolygon(pt(1),pt(2),Xelem(:,1),Xelem(:,2)))
            delInd(j) = true;
        end
    end
    channelLocNodes(:,delInd) = [];
    channelNodes(:,delInd) = [];
    channelNurbsParam(:,delInd) = [];
    channelNum(delInd) = [];
end
end

% this function finds the edge along which a line source may lie
% if the two end nodes of the line source are on the same edge, it returns the loc
% edge number. if not it returns 0
function channelLocEdge = channel_local_edge(parent)

locEdges(1,:) = parent.locEdge;
locEdges(2,1:3) = [3,1,2];
nChannels = numel(parent.channelNum);
channelLocEdge = zeros(nChannels,1);
for i = 1:nChannels
    locEdge = intersect(locEdges(:,parent.channelLocNodes(1,i)),...
                        locEdges(:,parent.channelLocNodes(2,i)),'stable');  
    locEdge(locEdge == 0) = [];                
    if (~isempty(locEdge))
        channelLocEdge(i) = locEdge;
    end
end
end

% this function removes a line source that intersect two ends of an edge
% but is not coincident with the edge
function parent = remove_non_coincident_channels(parent)
nChannels = numel(parent.channelNum);
del = false(nChannels,1);
for i = 1:nChannels
    if(parent.channelLocEdge(i) > 0 && parent.channelLocEdge(i) < 4)
        if(~parent.coincident(parent.channelLocEdge(i)))
            del(i) = true;
            %parent.isSpecial = false;
        end
    end
end


parent.channelNodes(:,del) = [];
parent.channelLocNodes(:,del) = [];
parent.channelNum(del) = [];
parent.channelNurbsParam(:,del) = [];
parent.channelLocEdge(del) = [];
if (numel(parent.channelNum) == 0)
    parent.type = 0;
end
end


% this function finds the collinear sets of nodes in a parent element
% (assumed to be triangular)
% the collinear sets are given by the terminal node local numbers instead of 
% the global numbers.  
%{
function collinearLoc = find_collinear_sets(parent,edge)
nEdge=3; % this is set to 3 because the parent element can only be a ...
         % triangular element
nNodes = numel(parent.nodes);
collinearLoc=zeros(nEdge,nNodes);
iEN=4:nNodes; % 4 because the parent element is assumed to be triangular
enrichNodeEdge=parent.edge(iEN);
for i=1:nEdge
    if (i==1)
        collinearLoc(i,1:2)=[1,2];
    elseif(i==2)
        collinearLoc(i,1:2)=[2,3];
    elseif(i==3)
        collinearLoc(i,1:2)=[3,1];
    end
    ind=(enrichNodeEdge==edge(i));
    nInd=nnz(ind);
    if(nInd>0)
        collinearLoc(i,3:2+nInd)=iEN(ind);
    end
end
% remove zero columns(s) of collinear
collinearLoc=collinearLoc(:,any(collinearLoc));
% remove row with <=2  non-zero entries
collinearLoc(sum(collinearLoc>0,2)<=2,:)=[];
end
%}
