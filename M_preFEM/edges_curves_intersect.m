%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Marcus Tan on 8/20/2013
%%% Last modified date: 1/16/2015
%%% Copyright 2013 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This function finds the intersection of the edges of the elements with
% nurbs curves. If multiple intersections per edge is found, local
% refinement by longest edge bisection is performed
% LIMITATION: Not yet able to handle branching point on channel edge
% INPUT:
% edge.edge_node
% node.coords
% node.n_node
% channels.nNurbs
% channels.nurbs
% channels.dim
% channels.kind
% tol.epsco
% tol.epsge
% tol.nurbsParam
% tol halfLineWidth: half the line segment width within which a terminal node of an edge is
% considered to intersect the line segment

% OUTPUT:
% edge.itrsect.num: number of intersections
% edge.itrsect.nSeg: a vector with itrsect.num entries specifying the number of nurbs
%               curve segments passing through a node 
% edge.itrsect.seg: a itrsect.num cell specifying the nurbs 
%               curve segments passing through a node 
% edge.itrsect.nurbsParam: a itrsect.num cell specfying the location on the
%                     nurbs curve segment that has the intersection point
% edge.itrsect.ipt: a vector with itrsect.num entries indicating the terminal
%                     point of a section of nurbs interface. if the
%                     intersection point is not a terminal point, it is set
%                      to zero.
% edge.itrsect.node: a vector with itrsect.num entries specifying the node number
% edge.itrsect.x
% edge.itrsect.y
% edge.itrsect.multItrsect: flag indicating whether an edge has multiple
%                           intersections and at least one intersection is 
%                           not an original node. if an edge belongs to an element
%                           with branching, this flag is false even if the
%                           edge has multiple intersections.
% node.n_node
% node.node_n
% node.coords

function [edge, node, elem, chanItrsectParams,count] ...
                    = edges_curves_intersect(edge, node, elem, ...
                                             channels, tol, refine, ...
                                             calcItrsectVel)
if (nargin < 6)
    refine.maxRefineLevel = 15;
    refine.refineJuncElem = false;
    %refine.maxAllowedAngleDiff = 30*pi/180;
end
%refine.cosMaxAllowedAngleDiff = cos(refine.maxAllowedAngleDiff);

fprintf('the max number of refinements is %i \n',refine.maxRefineLevel);
nEdge=size(edge.edge_node,1);
node.n_node = size(node.coords,1);
currentNode = node.n_node;

node.constraint.number = 0;
if(isfield(channels,'pt_temp'))
    ipt_temp = channels.pt_temp;
    node.constraint.temp_node = zeros(size(ipt_temp,1),1);
    node.constraint.temp_value = zeros(size(ipt_temp,1),1);
else
    ipt_temp = [];
    node.constraint.temp_node = [];
    node.constraint.temp_value = [];
end

% first stores the channel intersection parameters corresponding to the 
% enriched nodes. then they are combined with channel parameters 
% corresponding to original nodes
chanItrsectParams = cell(channels.nNurbs,1);
chanEnrichNodes = cell(channels.nNurbs,1);
chanOrigParams = cell(channels.nNurbs,1);
chanOrigNodes = cell(channels.nNurbs,1);

% compute derivative function of nurbs 
for j=1:channels.nNurbs
    channels.dnurbs(j) = nrbderiv(channels.nurbs(j));
    %channels.nurbsParams{j} = [];
    %channels.enrichNode{j} = [];
end
channels.dim = 2;


edge.itrsect = struct('num',0,'nSeg',cell(nEdge,1),'seg',cell(nEdge,1),...
                      'ipt',cell(nEdge,1),'nurbsParam',cell(nEdge,1),...
                      'node',cell(nEdge,1),'x',cell(nEdge,1),'y',cell(nEdge,1),...
                      'tangents',cell(nEdge,1),...
                      'potDirichlet',false,'multItrsect',false,'coincident',false);
if (calcItrsectVel)
    [edge.itrsect.vx,edge.itrsect.vy,edge.itrsect.designParamNum] = deal([]);
end

findItrsect = true;

% use a fast function to determine which edges have intersections when all
% the NURBS are linear
%
if (all(cat(1,channels.nurbs.order) == 2))
    checkEdgeIntersection ...
        = any(line_segment_intersect_indicator(...
                    cat(1,channels.lineSegs{:}),...
                    [node.coords(edge.edge_node(:,1),:), ...
                     node.coords(edge.edge_node(:,2),:)],...
                     tol.intersectEdges))';
else
    checkEdgeIntersection = true(nEdge,1); % indicates whether an edge has been modified or is new
end
%
%checkEdgeIntersection = true(nEdge,1);

count = 0; % refinement level

while (findItrsect)   
    % NOTE: When refinement happens, channels.nurbsParams and
    % channels.enrichNode only contain the information for edges whose
    % intersections have been recalculated
    for j=1:channels.nNurbs
        channels.nurbsParams{j} = [];
        channels.enrichNode{j} = [];
    end
    
    nElems = size(elem.elem_node,1);
    elemMarkers = false(nElems,1);
    elem.dualedge = elems_sharing_edge_nodes(elem.elem_node,node.n_node);    
    if refine.refineJuncElem
        % exclude branching points so that elements with branching points
        % are refined when multiple intersections occur
        BK4testing = channels.branch_kinks.XX(:,~channels.branch_kinks.isBranch);
        isBranch = false(1,size(BK4testing,2));
    else
       % include branching points so that elements with branching points 
        % in the interior are not refined even when multiple intersections
        % occur
        BK4testing = channels.branch_kinks.XX;
        isBranch = channels.branch_kinks.isBranch;
    end
   
    
    for i = 1:nEdge % for i  
        if(checkEdgeIntersection(i)) 
            edgeCoords = node.coords(edge.edge_node(i,:),:)';
            [edge.itrsect(i),currentNode,...
             chanItrsectParams,chanEnrichNodes, ...
             chanOrigParams,chanOrigNodes]...
                    = single_edge_curves_intersect(edge.itrsect(i), ...
                                                   edge.edge_node(i,:), ...
                                                   edgeCoords, currentNode, ...
                                                   channels, tol, ...
                                                   chanItrsectParams, ...
                                                   chanEnrichNodes, ...
                                                   chanOrigParams, ...
                                                   chanOrigNodes, ...
                                                   calcItrsectVel);                            
            checkEdgeIntersection(i) = false;    
        end

        % mark elements (excluding element with branching and 
        % element with a kink on its edge) sharing an edge
        % with multiple intersections per channel   
        seg = cell2mat(edge.itrsect(i).seg);
        
        %cosMaxAngleDiff = cosine_max_tangent_angle_differences(edge.itrsect(i).tangents);
        if (refine.maxRefineLevel ...
                && numel(seg) > numel(unique(seg)) ...
                &&  ~edge.itrsect(i).coincident)
            el1 = elem.dualedge(edge.edge_node(i,1),edge.edge_node(i,2));
            el2 = elem.dualedge(edge.edge_node(i,2),edge.edge_node(i,1));
            
            twoElemNodes = union(elem.elem_node(el1,:), ...
                                 elem.elem_node(el2,:),'stable');
 
            [~,sharedLocNodes] = intersect(elem.elem_node(el1,:),...
                                           elem.elem_node(el2,:),'stable');  
            twoElemCoords = node.coords(twoElemNodes,:);

            twoElemDT = delaunayTriangulation(twoElemCoords,sharedLocNodes');
            % determine whether there is a branching point in the interior
            % or a kink on the edge
            [elBKflags,BK4testing,isBranch] ...
                = interior_branching_edge_kink(BK4testing,...
                                               isBranch,...
                                               twoElemDT,...
                                               tol.halfLineWidth);   
            % do not refine elements sharing edge if one of the elements 
            % contain a branching point in its interior 
            % or a kink on its edge            
            if ~any(elBKflags)
                elemMarkers(el1) = true;
                elemMarkers(el2) = true;
                edge.itrsect(i).multItrsect = true;
            end          
        end
    end 
   
    % uncomment this block this see the mesh and intersections as the mesh
    % is refined
    %
    %plot_mesh_labels(node.coords,elem.elem_node,true,true,edge.edge_node)
    %{
    plot_mesh_curve_itrsect_junc(node,elem.elem_node,channels,...
                                 edge.itrsect,[],[],true,true,elemMarkers);
    box on
    set(gca,'XTick',[])
    set(gca,'YTick',[])
    set(gca,'XTickLabel',[])
    set(gca,'YTickLabel',[])
    %}
    if (count >= refine.maxRefineLevel)
        break
    end
    
    if (all(~elemMarkers))
        findItrsect = false;
    else
       % perform local edge refinemenet
       count = count+1;  
       fprintf('performing refinement number: %i\n',count)
       [node.coords,...
        elem.elem_node,...
        edge.edge_node,...
        checkEdgeIntersection,...
        elemProperty] ...
            = bisection(elemMarkers,...
                        node.coords,...
                        elem.elem_node,...
                        edge.edge_node,...
                        elem.dualedge,...
                        {elem.heatSource,elem.material});        
       elem.heatSource = elemProperty{1}';
       elem.material = elemProperty{2}';

       nEdge = size(edge.edge_node,1);
       
       nEdgeOld = numel(edge.itrsect); 
       % intersection information on modified edge is removed and that of
       % new edge is created
       for i = 1:nEdgeOld
            if(checkEdgeIntersection(i))
                edgeCoords = node.coords(edge.edge_node(i,:),:)';
                edgeLen = norm(edgeCoords(:,2)-edgeCoords(:,1),2);
                if (edgeLen < edge.minLength)
                    edge.minLength = edgeLen;
                end
                % remove the intersections in chanEnrichNodes, 
                % chanOrigNodes etc by searching through them.
                % not an efficient way but this will do for now
                for m = 1:numel(chanEnrichNodes)
                    [~,delInds] = intersect(chanEnrichNodes{m},edge.itrsect(i).node);
                    chanEnrichNodes{m}(delInds) = [];
                    chanItrsectParams{m}(delInds) = [];
                    [~,delInds] = intersect(chanOrigNodes{m},edge.itrsect(i).node);
                    chanOrigNodes{m}(delInds) = [];
                    chanOrigParams{m}(delInds) = [];
                end
                edge.itrsect(i) = clean_edge_itrsect(calcItrsectVel);             
            end
       end

       % initiate the new edges
       for i = (nEdgeOld+1):nEdge
            edge.itrsect(i) = clean_edge_itrsect(calcItrsectVel); 
       end
       
       %tol.halfLineWidth = edge.minLength*tol.halfLineWidthFrac;
       
       % old enrichment node labels before refinement

       oldEnNodes = cat(2,edge.itrsect(:).node);
       oldEnNodes = oldEnNodes(oldEnNodes > node.n_node);
       nOldEnNodes = numel(oldEnNodes);
       
       % create mapping from old enrichment node labels to new enrichment
       % node labels     
       node.n_node = size(node.coords,1);  
       if (nOldEnNodes > 0)
           %newEnNodes = sparse(oldEnNodes,1,...
           %                    node.n_node+ (1:nOldEnNodes),...
           %                    max(oldEnNodes),1);
           newEnNodes = 1:max(oldEnNodes);
           newEnNodes(oldEnNodes) = node.n_node+ (1:nOldEnNodes);
           for i = 1:nEdgeOld
               if edge.itrsect(i).num
                    edge.itrsect(i).node = newEnNodes(edge.itrsect(i).node);
               end              
           end
           for m = 1:numel(chanEnrichNodes)
                chanEnrichNodes{m} = newEnNodes(chanEnrichNodes{m});
           end
           currentNode = node.n_node + nOldEnNodes;
       end
    end
    
   
    %   
end


node.nOriginalNode = node.n_node;


for i=1:nEdge
    for j = 1:edge.itrsect(i).num
        % add enrichment nodes
        if (edge.itrsect(i).node(j) > node.n_node)
             node.coords(edge.itrsect(i).node(j),1) = edge.itrsect(i).x(j);
             node.coords(edge.itrsect(i).node(j),2) = edge.itrsect(i).y(j);
        end
        % check if any of the intersection points are the points where the
        % temperature are specified. If so update the intersection point as a
        % Diricihlet node
        if(~isempty(ipt_temp) && edge.itrsect(i).potDirichlet)              
            iptInd = ipt_temp(:,1)==edge.itrsect(i).ipt(j);
            if(any(iptInd))

                %[node.Dirichlet,addedDirichlet] = add_Dirichlet_node(node.Dirichlet,...
                %                                                     edge.itrsect(i).node(j),...
                %                                                     ipt_temp(iptInd,2));                         
                node.constraint.number = node.constraint.number+1;
                fprintf('added constraint node %i\n', edge.itrsect(i).node(j));
                node.constraint.temp_node(node.constraint.number)...
                    = edge.itrsect(i).node(j);
                fprintf('constraint value = %g\n',ipt_temp(iptInd,2));
                node.constraint.temp_value(node.constraint.number)... 
                    = ipt_temp(iptInd,2);
                ipt_temp(iptInd,:) = [];                
            end  
        end
    end    
end %i

edge.itrsect = rmfield(edge.itrsect,{'potDirichlet'});
if(~isempty(ipt_temp))
    warning('Unable to assign some temperature points on channel to any node, reduce nurbs tolerance')
    disp(ipt_temp)
end
node.n_node = size(node.coords,1);

elem.n_elem = size(elem.elem_node,1);
                    
% sort channel intersection parameters so that elem_branching_n_kinks.m
% can search through them faster
for j=1:channels.nNurbs
    chanItrsectParams{j} = union(chanItrsectParams{j},chanOrigParams{j});
end


end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Marcus Tan on 6/27/2014
%%% Last modified date: 6/27/2014
%%% Copyright 2014 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function itrsect = clean_edge_itrsect(calcItrsectVel)
    itrsect.num = 0;
    itrsect.nSeg = [];
    itrsect.seg = [];
    itrsect.ipt = [];
    itrsect.nurbsParam = [];
    itrsect.node = [];
    itrsect.x = [];
    itrsect.y = [];
    itrsect.tangents = [];
    itrsect.potDirichlet = false;
    itrsect.multItrsect = false;
    itrsect.coincident = false;
    if (calcItrsectVel)
        itrsect.vx = [];
        itrsect.vy = [];
        itrsect.designParamNum = [];
    end
    
end
function [elBKflags,BK4testing,isBranch] ...
                = interior_branching_edge_kink(BK4testing,...
                                               isBranch,...
                                               DT,...
                                               halfLineWidth)
elBKflags = false(2,1);    

[elemNums,baryCoords] = pointLocation(DT,BK4testing');
delInds = false(1,size(BK4testing,2));
for i = 1:numel(elemNums)
   if isnan(elemNums(i))
       continue
   else
       %locEdge = loc_edge_from_bary_coords(baryCoords(i,:),tolBary);
       elemCoords = DT.Points(DT.ConnectivityList(elemNums(i),:),:);
       [in,locEdge] = inTriangle(BK4testing(1,i),BK4testing(2,i),...
                                 elemCoords(:,1),elemCoords(:,2),halfLineWidth); 
       if isBranch(i) && locEdge == 0
           % branching point in element interior
           elBKflags(elemNums(i)) = true;
           % remove the branching point from further testing
           delInds(i) = true;
       elseif ~isBranch(i) && locEdge > 0
           % kink on element edge
           elBKflags(elemNums(i)) = true;
       end
   end
end
BK4testing(:,delInds) = [];
isBranch(delInds) = [];
end

%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Marcus Tan on 7/2/2014
%%% Last modified date: 7/13/2014
%%% Copyright 2014 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check whether any junction in junc4test is in an element. return juncFlag
% = true, if it is.  
% Also, it is removed from the list junc4test and the list is returned.
function [juncFlag,junc4test] = junc_in_single_element(juncFlag,...
                                                       junc4test,...
                                                       pts,...
                                                       elem_coords,...
                                                       halfLineWidth)
if (~juncFlag)          
    for i = 1:numel(junc4test)
        j = junc4test(i);
        %if (inpolygon(pts(j,1),pts(j,2),elem_coords(:,1),elem_coords(:,2)))
        if (inTriangle(pts(j,1),pts(j,2),elem_coords(:,1),elem_coords(:,2),halfLineWidth))
            junc4test(i) = 0;
            juncFlag = true;
        end
    end
    junc4test = junc4test(junc4test > 0);
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Marcus Tan on 7/13/2014
%%% Last modified date: 11/7/2014
%%% Copyright 2014 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check whether any kink in kink_coords is in an element. return kinkFlag
% = true, if it is.  
function kinkFlag = kink_in_single_element(kinkFlag,...
                                           kink_coords,...
                                           elem_coords,...
                                           halfLineWidth)
if (~kinkFlag)
    nKinks = size(kink_coords,2);
    for i = 1:nKinks
        if (inTriangle(kink_coords(1,i),kink_coords(2,i),elem_coords(:,1),elem_coords(:,2),halfLineWidth))
            kinkFlag = true;
        end
    end
end
end
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Marcus Tan on 7/13/2014
%%% Last modified date: 7/13/2014
%%% Copyright 2014 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check whether any kink in kink_coords is in an element. return kinkFlag
% = true, if it is.  
% Also, it is removed from the list kink_coords and the list is returned.
%{
function [kinkFlag,kink_coords] = kink_on_single_element_edge(kinkFlag,...
                                                              kink_coords,...
                                                              elem_coords,...
                                                              halfLineWidth)
if (~kinkFlag)
    nKinks = size(kink_coords,2);
    mask = false(nKinks,1);
    for i = 1:nKinks
        [~,locEdge] = inTriangle(kink_coords(1,i),kink_coords(2,i),elem_coords(:,1),elem_coords(:,2),halfLineWidth); 
        if (locEdge)
            mask(i) = true;
            kinkFlag = true;
        end
    end
    kink_coords(:,mask) = [];
end
end
%}
%{
function flag = itrsect_near_parallel2edge(itrsect,edgeCoords,channels,cosAngleTol)
if (itrsect.num ~= 1)
    flag = false;
    return
end
flag = false;
edgeVec = edgeCoords(:,2) - edgeCoords(:,1);
if (itrsect.nSeg == 1)
    j = itrsect.seg{1};
    [~,tangent] = nrbdeval(channels.nurbs(j),channels.dnurbs(j),itrsect.nurbsParam{1});
    cosAngle = abs(dot(edgeVec,tangent(1:2))/(norm(edgeVec)*norm(tangent(1:2))));
    if (cosAngle > cosAngleTol)
        flag = true;
        return;
    end
end
end
%}
%%% Created by Marcus Tan on 1/5/2014
%%% Last modified date: 3/31/2014
%%% Copyright 2014 University of Illinois at Urbana-Champaign. All rights reserved
% this function calculates and sort a set of points (xs(1),ys(1)),...,(xs(n),ys(n))
% using the distance from a specified point (xo,yo)
%{
function [minInd,minDist] = sort_pts_based_on_distance(xs,ys,xo,yo)
minDist = inf;
minInd = 0;
if(isempty(xs) || isempty(ys))
    return
end
npts = numel(xs);

for i = 1:npts
    dist = norm([xs(i)-xo, ys(i)-yo],2);
    if (dist < minDist)
        minDist = dist;
        minInd = i;
    end
end
end
%}

%{
% this function determines the absolute value of cosine of the 
% maximum difference of the angles between the tangents at the intersections
function cosMaxAngleDiff = cosine_max_tangent_angle_differences(tangents)
nTangents = size(tangents,2);
cosMaxAngleDiff = inf;
for i = 1:nTangents-1;
    for j = i+1:nTangents
        cosAngle = abs(dot(tangents(:,i),tangents(:,j))/(norm(tangents(:,i))*norm(tangents(:,j))));
        if (cosAngle < cosMaxAngleDiff)
            cosMaxAngleDiff = cosAngle;
        end
    end
end
end
%}