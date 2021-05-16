%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Marcus Tan on 5/27/2013
%%% Last modified date: 10/28/2014
%%% Copyright 2013 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This function creates child elements from a parent element
% REMARKS:
%   i) When there is a branching point, a special element where the same
%       channel coincide with two edges of the element and multiple
%       channels but not branching point in an element, triangular child
%       elements are not merged to form quadrilateral child elements
%   ii) When there is a branching point, the triangular child elements are
%       such that each child element cannot have more than one channel 
%   iii) Triangular elements are not merged when there is multiple channels
%        in a parent element or when polynomial IGFEM is used 
% INPUT:
% parent:
% elem_material:
% material:
% node:
% slenderTol:
% 
% OUTPUT:
% parent(i).child(j).isTriangle
% parent(i).child(j).nodes: the global nodes of a child element.
%                           the exterior nodes are arranged in CCW order
%                           first followed by the interior nodes
% parent(i).child(j).locPaNodes: column vector of child local node number wrt parent
%                                nodes. NOTE: DOES NOT INCLUDE INTERIOR
%                                ENRICHMENT NODES
% parent(i).child(j).locPaEnNodes:column vector of child enrichment node
%                                 wrt parent nodes. NOTE: DOES NOT INCLUDE INTERIOR
%                                ENRICHMENT NODES
% parent(i).child(j).locEnNodes: column vector of child enrichment node wrt
%                                 its own nodes. NOTE: DOES NOT INCLUDE INTERIOR
%                                ENRICHMENT NODES
% parent(i).child(j).nurbsSeg: vector of nurbs data structure 
%                              describing the nurbs curve segments of the
%                              channels
% parent(i).child(j).dnurbsSeg: derivative of nurbsSeg  
% parent(i).child(j).channelNum: column vector of channel labels. must be
%                                of the same length as the number of columns 
%                                in channelNodes
% parent(i).child(j).channelNodes:  nChannel cell array specifying the 
%                                   enrichment nodes of the channels
%                                   for ex. channelNodes = [{1,2},{4,5,6}] 
%                                   means that two edges of the child with 
%                                   nodes 1,2 and 4,5,6 are channels.
%                                   The enrichment nodes are ordered in 
%                                   direction of increasing NURBS parameters
% parent(i).child(j).channelLocNodes: refers the local node number (REMOVED)
% parent(i).child(j).channelNurbsParam: 2 x nChannels matrix of the nurbs 
%                                       parametric coordinates of the channel 
%                                       terminal nodes
% parent(i).child(j).nurbsSurf: nurbs data structure
%                               describing the surface of the child elements 
% parent(i).child(j).uv: vector of the fixed parameter direction
%                        of the channels. 
%                        1: u-direction in parametric space, 
%                        2: v-direction in parametric space
% parent(i).child(j).uvParam: the value of the fixed parametric coordinate in the
%                             child.uv-direction
%                             Ex 1. nurbs segment runs along u direction at v=0, then
%                                   uv=2, uvParam=0
%                             Ex 2. nurbs segment runs along v direction at u=1, then
%                                   uv=1, uvParam=1
% parent(i).child(j).nurbsInd: cell arrays of the indices of the
%                              nurbs.coef(3,:,:) or the nurbs basis-function of
%                              corresponding to the nodes of a child element 
%                              Ex. child.node = [3,5,6]
%                              nodes 5 and 6 are intersection nodes
%                              nurbsInd = {[],[1;1],[2,1;2,2]'} means that 
%                              the nurbs basis functions used for enrichment for nodes
%                              5 and 6 are N11,N21+N22
%                       
function [parent,node]=child_elements_nurbs(parent,...
                                            node, ...
                                            collinearTol, ...
                                            slenderTol, ...                                            
                                            polyIGFEM)
if (polyIGFEM)
    mergeTriangles = false;
    fprintf('\nNOTE: triangles are not merged in child_elements_nurbs\n')
else
    mergeTriangles = true;
end

if (isstruct(slenderTol))
    cosMinAngleTol = cos(slenderTol.minAngle*pi/180.0); 
end
nNodesPerElem = 3; % no. of nodes per original element

childFields2remove = {'map2parentChannel'};

nElems = length(parent);
localNodeSets = [1,2;3,1;2,3]; % for edge_flipping in the case where there is branching and child element has 2 line sources
for p = 1:nElems
   
    if(parent(p).type > 1)    
        pts = node.coords(parent(p).nodes,:);
        DT = delaunayTriangulation(pts,parent(p).channelLocNodes');
        TR = DT.ConnectivityList;
        nTri=size(TR,1);
        delTR = false(nTri,1);
        % remove the points from triangulation when they are supposed to be
        % collinear
        %{
        nColSets = size(parent(p).collinearLoc,1);
        for i = 1:nTri
            for k = 1:nColSets
                diff=setdiff(TR(i,:) ,parent(p).collinearLoc(k,:));
                if(isempty(diff))
                    delTR(i) = true;
                end
            end
        end
        %}
        % remove the points from triangulation if they are collinear to
        % within a tolerance
        for i = 1:nTri
            xv = pts(TR(i,:),1);
            yv = pts(TR(i,:),2);
            if polyarea(xv,yv) < collinearTol
                delTR(i) = true;
            end
        end
        TR(delTR,:)=[];
 
        nChannels = numel(parent(p).channelNum);
        rec = [];
        if (nChannels < 2 && mergeTriangles)
            % merge 2 triangles to form a quadrilatel if 
            % (i) the triangles do not share a constraint edge.
            % (ii) the resulting quadrilateral do not have an interior angle> 180 deg.
            [TR,rec] = merge_triangles(TR,pts,parent(p).channelLocNodes);       
        end

        locNodes = [[TR,zeros(size(TR,1),1)];rec];
        %locNode=[padarray(TR,[0 1],0,'post');Rec];
        
        nChildren=size(locNodes,1);
        approxnurbs = false;
        
        parent(p).child = struct('isTriangle',false,...
                                'nodes',cell(nChildren,1), ...
                                'locPaNodes',cell(nChildren,1), ...
                                'locPaEnNodes',cell(nChildren,1), ...
                                'locEnNodes',cell(nChildren,1), ...
                                'channelNum',cell(nChildren,1), ...
                                'channelNodes',cell(nChildren,1), ...
                                'channelLocNodes',cell(nChildren,1), ...
                                'channelNurbsParam',cell(nChildren,1), ...
                                'nurbsSeg',cell(nChildren,1), ...
                                'dnurbsSeg',cell(nChildren,1), ...
                                'nurbsSurf',cell(nChildren,1),...
                                'dnurbsSurf',cell(nChildren,1),...
                                'uv',cell(nChildren,1),...
                                'uvParam',cell(nChildren,1),...
                                'nurbsInd',cell(nChildren,1),...
                                'hullx',cell(nChildren,1),...
                                'hully',cell(nChildren,1));
        for i=1:nChildren
            nonZeros = locNodes(i,:)>0;
            parent(p).child(i).locPaNodes = locNodes(i,nonZeros)';
            parent(p).child(i).nodes = parent(p).nodes(parent(p).child(i).locPaNodes);
            parent(p).child(i).locEnNodes ...
                = find(parent(p).child(i).locPaNodes > nNodesPerElem); 
            parent(p).child(i).locPaEnNodes ...
                = parent(p).child(i).locPaNodes(parent(p).child(i).locEnNodes);
            %parent(p).child(i).conductivity = material(elem_material(p)).conductivity;
            if(nnz(nonZeros)>3)
                parent(p).child(i).isTriangle=false;
            else
                parent(p).child(i).isTriangle=true;
            end
            if (~polyIGFEM)
                % check for slender element when the check has not be done, the
                % tolerance is specified and when any nurbsSeg is order 3 and above
                if(~approxnurbs && isstruct(slenderTol) && any(cat(2,parent(p).nurbsSeg.order) > 2))
                    if(parent(p).child(i).isTriangle)
                        cosMinAngle = cosine_min_angle_2D_triangle(node.coords(parent(p).child(i).nodes,:)'); 
                         if(cosMinAngle > cosMinAngleTol)
                            disp(['child ',num2str(i),' of parent element ',num2str(p),' has min angle ',num2str(acos(cosMinAngle)*180/pi)])
                            approxnurbs = true;
                         end
                    else
                         ratio = aspect_ratio_2D(node.coords(parent(p).child(i).nodes,:)');
                         if (ratio >  slenderTol.maxAspectRatio)
                            disp(['child ',num2str(i),' of parent element ',num2str(p),' has aspect ratio ',num2str(ratio)])
                            approxnurbs = true;
                         end
                    end               
                end
            end
        end      
        

        % if any child element has an aspect ratio greater than
        % maxAspectRatio, approximate NURBS curve with a line segment
        if approxnurbs
            disp('replacing NURBS curves with line segments')
            for i=1:numel(parent(p).channelNum)
                parent(p).nurbsSeg(i) = nurbs_curve_fit_segment(parent(p).nurbsSeg(i),[0,1],2,2);
                parent(p).dnurbsSeg(i) = nrbderiv(parent(p).nurbsSeg(i));
            end
        end
        % find channels of each child element
        parent(p) = child_elements_channels(parent(p),polyIGFEM);
        
        % if element has a branching point and a child element has > 1 line source,  
        % swap edges until only there is only at most 1 line source per
        % child element
        for c = 1:numel(parent(p).child)
            % swap the edge of ajacent child elements if one child
            % element has two line sources
            nChannels = numel(parent(p).child(c).channelNum); 
            if (nChannels == 2)
                % lsEdge = 3 if local nodes are 1 and 2, 5 if local
                % nodes are 2 and 3 and 4 if local nodes are 1 and 3.
                % determines which edge is to be swapped
                lsEdge1 = parent(p).child(c).channelLocNodes(1,1)... 
                         +parent(p).child(c).channelLocNodes(2,1);
                lsEdge2 = parent(p).child(c).channelLocNodes(1,2)... 
                         +parent(p).child(c).channelLocNodes(2,2);
                nonLsEdge = setdiff([3,5,4],[lsEdge1,lsEdge2],'stable');
                % nodes of the edge to be swapped
                edgeLocNodes = parent(p).child(c).locPaNodes(localNodeSets(nonLsEdge-2,:))';
                % opposite child element sharing the edge to be swapped
                % can make more efficient later !!!!!!!!!!!
                [row1,~] = find(locNodes == edgeLocNodes(1));
                [row2,~] = find(locNodes == edgeLocNodes(2));
                c2 = setdiff(intersect(row1,row2),c,'stable');
                % check that quadrilateral formed by combining the the two
                % triangular elements is convex before performing edge
                % swapping
                if(~isempty(c2))
                    quadLocNodes = union(parent(p).child(c).locPaNodes,...
                                         parent(p).child(c2).locPaNodes,...
                                         'stable'); 
                                   
                    quadPts = pts(quadLocNodes,:); 

                    hullInd = convhull(quadPts(:,1),quadPts(:,2)); 
                    
                    if (numel(hullInd) < 5)
                        swapEdge = false;
                    else
                        swapEdge = true;
                    end
                else
                    swapEdge = false;
                end
                if(swapEdge)
                    adjChildrenLocNodes = [parent(p).child(c).locPaNodes(1:3)'; ...
                                            parent(p).child(c2).locPaNodes(1:3)'];
                    [~,newChildrenLocNodes] = flip_edge(edgeLocNodes,adjChildrenLocNodes);                    
                    parent(p).child(c).locPaNodes = newChildrenLocNodes(1,:)';
                    parent(p).child(c2).locPaNodes = newChildrenLocNodes(2,:)';
                    locNodes(c,1:3) = parent(p).child(c).locPaNodes;
                    locNodes(c2,1:3) = parent(p).child(c2).locPaNodes;       
                    parent(p).child(c).nodes = parent(p).nodes(newChildrenLocNodes(1,:));
                    parent(p).child(c2).nodes = parent(p).nodes(newChildrenLocNodes(2,:));
                    parent(p).child(c).locEnNodes ...
                        = find(parent(p).child(c).locPaNodes > nNodesPerElem); 
                    parent(p).child(c).locPaEnNodes ...
                        = parent(p).child(c).locPaNodes(parent(p).child(c).locEnNodes);
                    parent(p).child(c2).locEnNodes ...
                        = find(parent(p).child(c2).locPaNodes > nNodesPerElem); 
                    parent(p).child(c2).locPaEnNodes ...
                        = parent(p).child(c2).locPaNodes(parent(p).child(c2).locEnNodes);
                    
                    % can make more efficient by updating the channel info
                    % of only the two child elements that have been
                    % affected
                    parent(p) = child_elements_channels(parent(p),polyIGFEM);        
                else
                    if (parent(p).child(c).isTriangle)
                        fprintf('warning: 2 line sources in element %i, triangular child %i \n',p,c)
                    else
                        fprintf('warning: 2 line sources in element %i, quadrilateral child %i \n',p,c)
                    end
                end
            elseif (nChannels > 2)
                fprintf('warning: > 2 line sources per child element in element %i\n',p)
            end
        end

        
        %parent(p).child = rmfield(parent(p).child,'locPaNodes'); 
        if (~polyIGFEM)
            for i=1:numel(parent(p).child) 

                % the child element does not have any line source
                if (numel(parent(p).child(i).channelNum) == 0)
                    % this child element can have at most one enrichment node
                    % put this enrichment node as the first node
                    ind = find(parent(p).child(i).nodes > node.nOriginalNode,1,'first');
                    if (~isempty(ind))
                        rshift = numel(parent(p).child(i).nodes)-ind+1;
                        lshift = 1-ind;
                        if (rshift > abs(lshift))
                            rshift = lshift;
                        end
                        parent(p).child(i).nodes = circshift(parent(p).child(i).nodes,[0,rshift]);
                    end
                    xel=node.coords(parent(p).child(i).nodes,:)';
                    parent(p).child(i).channelNodes = [];
                    parent(p).child(i).channelLocNodes = [];
                    parent(p).child(i).channelNum = [];
                    parent(p).child(i).nurbsSeg = [];
                    parent(p).child(i).dnurbsSeg = [];
                    [parent(p).child(i).nurbsSurf,~,~]...
                               =child_elem_nurbs_surface(xel,...
                                                         parent(p).child(i).channelLocNodes,...
                                                         parent(p).child(i).nurbsSeg);
                % the child element have at least one line source
                else
                    % ASSUMPTION: at most two line source or nurbs segments per child element
                    xel=node.coords(parent(p).child(i).nodes,:)';
                    [parent(p).child(i).nurbsSurf,...
                     parent(p).child(i).uv,...
                     parent(p).child(i).uvParam]...
                               =child_elem_nurbs_surface(xel,...
                                                         parent(p).child(i).channelLocNodes,...
                                                         parent(p).child(i).nurbsSeg); 

                end
                parent(p).child(i).dnurbsSurf=nrbderiv(parent(p).child(i).nurbsSurf);

                % convex hull of a child element for interpolate_solution            
                xChild = node.coords(parent(p).child(i).nodes,1);
                yChild = node.coords(parent(p).child(i).nodes,2);
                hullInd = convhull(xChild,yChild);
                parent(p).child(i).hullx = xChild(hullInd);
                parent(p).child(i).hully = yChild(hullInd);

            end        
              
            [parent(p),new_coords] ...
                    = create_interior_enrichment_nodes(parent(p),...
                                                       size(node.coords,1));   
            % update locPaNodes, locEnNodes, locPaEnNodes

            if (~isempty(new_coords))           
                for i = 1:numel(parent(p).child)
                    [~,parent(p).child(i).locPaNodes] ...
                        = ismember(parent(p).child(i).nodes, ...
                                   parent(p).nodes);
                    parent(p).child(i).locEnNodes ...
                        = find(parent(p).child(i).locPaNodes > nNodesPerElem); 
                    parent(p).child(i).locPaEnNodes ...
                        = parent(p).child(i).locPaNodes(parent(p).child(i).locEnNodes);
                end
                node.coords = [node.coords;new_coords];
                node.n_node = size(node.coords,1);
            end

            parent(p) = nurbs_shape_func_ind(parent(p),node.nOriginalNode);
        end
        parent(p).child = rmfield(parent(p).child,childFields2remove);      
    end
end

end

% this function finds the channels of the child elements in a parent
% element
% INPUT:
%   1 <= i <= nChildren
%   parent.channelNodes: 2 x number of channels matrix
%   parent.child(i).nodes: column vector of child global node number
% OUTPUT:
%   1 <= i <= nChildren
%   parent.child(i).channelNum
%   parent.child(i).channelNodes
%   parent.child(i).channelLocNodes (REMOVED)
%   parent.child(i).nurbsSeg
%   parent.child(i).dnurbsSeg
%   parent.child(i).channelNurbsParam
%   parent.child(i).map2parentChannel 
function [parent,child2parentChannel] = child_elements_channels(parent,polyIGFEM)
nChildren = numel(parent.child);
nChannels = numel(parent.channelNum);
%nSharedChildren = zeros(nChannels,1);
child2parentChannel = zeros(nChannels,nChildren);
for i = 1:nChildren
    m = 0;
    parent.child(i).channelNodes = cell(nChannels,1);
    parent.child(i).channelLocNodes = zeros(2,nChannels);
    parent.child(i).channelNurbsParam = zeros(2,nChannels);
    parent.child(i).channelNum = zeros(nChannels,1);
    parent.child(i).map2parentChannel = zeros(nChannels,1); 
    for j=1:numel(parent.channelNum)
        if(isempty(setdiff(parent.channelNodes(:,j),parent.child(i).nodes)))
            m = m + 1;
            parent.child(i).map2parentChannel(m) = j;
            parent.child(i).channelNodes{m}...
                                = parent.channelNodes(:,j);  
            
            parent.child(i).channelLocNodes(:,m) ...
                        = [find(parent.child(i).nodes == parent.channelNodes(1,j),1,'first'),...
                           find(parent.child(i).nodes == parent.channelNodes(2,j),1,'first')];

            parent.child(i).channelNurbsParam(1:2,m) = parent.channelNurbsParam(1:2,j); 
            parent.child(i).channelNum(m) = parent.channelNum(j);
            if (~polyIGFEM)
                parent.child(i).nurbsSeg(m) = parent.nurbsSeg(j);
                parent.child(i).dnurbsSeg(m) = parent.dnurbsSeg(j);
            end
        end
    end
    parent.child(i).channelNodes = parent.child(i).channelNodes(1:m);
    parent.child(i).channelLocNodes = parent.child(i).channelLocNodes(1:2,1:m);
    parent.child(i).channelNurbsParam = parent.child(i).channelNurbsParam(1:2,1:m);
    parent.child(i).channelNum = parent.child(i).channelNum(1:m);
    parent.child(i).map2parentChannel = parent.child(i).map2parentChannel(1:m); 
    %parent.child(i).nurbsSeg(m:end) = [];
    %parent.child(i).dnurbsSeg(m:end) = [];
end

   
end

% this function add enrichment nodes corresponding to the interior control
% points of the nurbs line source segments 
% Remarks: if interior enrichment nodes are added at an edge of the parent
% element, there is no way to prevent the same enrichment nodes from being
% added in the adjacent child element in the adjacent parent element. But
% this cannot happen unless the NURBS segment at the edge is of order
% greater than 2. However, when the NURBS segment is of order greater than
% 2, it is very unlikely that it will coincide with the edge. 
function [parent, ...
          new_coords] ...
                = create_interior_enrichment_nodes(parent,...
                                                   nNodes)
    nChannels = numel(parent.channelNum);
    nChildren = numel(parent.child);
    new_coords = zeros(nChannels*max(cat(1,parent.nurbsSeg.number)-2),2);
    nNewCoords = 0;
    child_nodes = zeros(nChildren,4);                                           
    for i = 1:nChildren  
        if (parent.child(i).isTriangle)
            child_nodes(i,1:3) = parent.child(i).nodes;
        else
            child_nodes(i,:) = parent.child(i).nodes;
        end
    end
    for i = 1:nChannels
        if (parent.nurbsSeg(i).number > 2)
            % find childs sharing a line source
            [childInds,~] = find(child_nodes == parent.channelNodes(1,i));
            [childInds2,~] = find(child_nodes == parent.channelNodes(2,i));
            childInds = intersect(childInds,childInds2,'stable');
            % coordinates of new interior enrichment nodes
            for j = 2:parent.nurbsSeg(i).number-1
                nNewCoords = nNewCoords + 1;
                new_coords(nNewCoords,:) ...
                        = parent.nurbsSeg(i).coefs(1:2,j)'./parent.nurbsSeg(i).coefs(4,j); 
            end
            nIntEnNodes = parent.nurbsSeg(i).number-2;
            for j = 1:length(childInds);
                parent.child(childInds(j)).nodes(end+1:end+nIntEnNodes) ...
                        = (nNodes+1):(nNodes+nIntEnNodes);
                indLS = (parent.child(childInds(j)).map2parentChannel == i);
                parent.child(childInds(j)).channelNodes{indLS}(end+1:end+nIntEnNodes) ...
                        = (nNodes+1):(nNodes+nIntEnNodes);
            end
            parent.nodes(end+1:end+nIntEnNodes) = (nNodes+1):(nNodes+nIntEnNodes);
            nNodes = nNodes + nIntEnNodes;
        end        
    end
    
    % place the 2nd terminal node of a channel at the end
    for i = 1:nChildren   
        for j = 1:numel(parent.child(i).channelNum)
            [parent.child(i).channelNodes{j}(end),parent.child(i).channelNodes{j}(2)] ...
                    = deal(parent.child(i).channelNodes{j}(2),parent.child(i).channelNodes{j}(end));
        end
    end
    
    new_coords = new_coords(1:nNewCoords,:);
end

% this function defines the nurbs shape functions indices in the nurbs surface
% description of a child element
function parent = nurbs_shape_func_ind(parent,nOriginalNodes)
for i = 1:numel(parent.child)
    parent.child(i).nurbsInd = cell(1,numel(parent.child(i).nodes));
    if (parent.child(i).isTriangle)
        maxLocNode = 3;
    else
        maxLocNode = 4;
    end
    
    for j = 1:numel(parent.child(i).channelNum)
        [~,indInChildNodes] = ismember(parent.child(i).channelNodes{j},...
                                       parent.child(i).nodes);                        
        isCCWflag = isCCW(parent.child(i).channelLocNodes(:,j),maxLocNode);
        nLSNodes = length(parent.child(i).channelNodes{j}); % must be same length as indInChildNodes
        % if the first line source is not CCW, its nurbsSeg is reversed
        % when forming the NURBS surface
        % if the second line source is not CCW, its nurbsSeg is NOT
        % reversed when forming the NURBS surface
        if ((j == 1 && isCCWflag) || (j == 2 && ~isCCWflag))
            for k = 1:nLSNodes
                  % ind = (vInd-1)*nurbsSurf.number(1)+uInd; 
                  parent.child(i).nurbsInd{indInChildNodes(k)}(end+1) ...
                      = sub2ind(parent.child(i).nurbsSurf.number,k,j);
            end
        else
            for k = 1:nLSNodes
                  % ind = (vInd-1)*nurbsSurf.number(1)+uInd;
                  parent.child(i).nurbsInd{indInChildNodes(k)}(end+1) ...
                      = sub2ind(parent.child(i).nurbsSurf.number,nLSNodes+1-k,j);
            end
        end
    end
    
    if (numel(parent.child(i).channelNum) > 0)
        % remove nurbs shape functions for original nodes
        for j = 1:numel(parent.child(i).nodes)
            if (parent.child(i).nodes(j) <= nOriginalNodes)
                parent.child(i).nurbsInd{j} = [];
            end
        end    
    else
        % even when there is no line source, enrichment function has to be
        % introduced if there is an enrichment node. Note that it's
        % impossible to have more than one enrichment node and that I have
        % assumed that if an enrichment node is present, it is the first
        % node. 
        if (parent.child(i).nodes(1) > nOriginalNodes)
            parent.child(i).nurbsInd{1} = [1,1];
        end
    end
end
end

