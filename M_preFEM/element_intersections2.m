%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Marcus Tan on 7/3/2014
%%% Last modified date: 7/23/2014
%%% Copyright 2014 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find the intersections of an element with the interface/channel based on
% the intersections between the interface an the edges
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
%                    possible.  (REMOVED)
% parent(i).isSpecial: an element is a special element if
%                   a single line source coincide with 2 edges of an
%                   element  (REMOVED)
% parent(i).multItrsect: true if there is more than 2 intersections per edge
%                     and at least one intersection is not an original
%                     node. if an edge belongs to an element with branching, 
%                     this flag is false even if the edge has multiple intersections.
% parent(i).nItrsect: number of intersection points
% parent(i).nodes: 1D array of the global node numbers (both enrichment nodes 
%                  and original nodes) of the parent element. 
%                  The nodes on the edges are arranged in counter-clockwise order
%                  Interior nodes are listed last in no particular order.
% parent(i).edge: a length parent.nNode array of edges upon which the nodes
%                 lie. For nodes at the vertices, the convention of the edge
%                 between vertex 1 and 2 is assigned to vertex 1, the edge
%                 between vertex 2 and 3 being assigned to vertex 2 and the
%                 edge between vertex 3 and 1 being assigned to vertex 3 is
%                 applied.
% parent(i).nSeg: a length parent.nNode array of the number of nurbs curve 
%                 segments passing through the nodes of the corresponding
%                 indices. nSeg is 0 if a node is not an intersection point.
% parent(i).seg: a length parent.nNode array of cells of the nurbs segments
%                passing through the nodes of the corresponding indices. seg
%                is an empty cell if a node is not an intersection point
% parent(i).nurbsParam: a length parent.nNode array of cells of the nurbs parametric 
%                    coordinate of the intersection point.
% parent(i).coincident: boolean vector indicating if an edge has a
%                       coindicent channel
function parent = element_intersections(elem_nodes,...
                                        dualedge,...
                                        edge_nodes,...
                                        itrsect)
nElems = size(elem_nodes,1);
parent = struct('type',0, ...
                'multItrsect',false, ...
                'nItrsect',0, ...
                'coincident',false(1,3),...
                'nodes',cell(nElems,1), ...
                'edge',cell(nElems,1),...
                'locEdge',cell(nElems,1),...
                'nSeg',cell(nElems,1), ...
                'seg',cell(nElems,1),...
                'nurbsParam',cell(nElems,1),...
                'conductivity',cell(nElems,1));         
for i=1:size(edge_nodes,1)
    if(itrsect(i).num>0)
        % find the elements sharing the edge given by ii
        adjElem = [dualedge(edge_nodes(i,1),edge_nodes(i,2)),...
                   dualedge(edge_nodes(i,2),edge_nodes(i,1))];
        adjElem = adjElem(adjElem > 0);
        % update candidate parent elements with the intersection nodes and the
        % locations  (in the parametric domain) of the nurbs curve intersected
        % by the element edges
        for j = 1:length(adjElem)
            %if (any(ismember(itrsect(i).ipt,channelJunc)))
            %    parent(adjElem(j)).hasJunc = true;
            %end
            k1 = parent(adjElem(j)).nItrsect + 1;            
            k2 = parent(adjElem(j)).nItrsect + itrsect(i).num;
            parent(adjElem(j)).nItrsect = k2;
            parent(adjElem(j)).nodes(k1:k2) = itrsect(i).node;
            parent(adjElem(j)).edge(k1:k2) = repmat(i,1,itrsect(i).num);               
            parent(adjElem(j)).nSeg(k1:k2) = itrsect(i).nSeg;
            parent(adjElem(j)).seg(k1:k2) = itrsect(i).seg;

            locEdge = local_edge(elem_nodes(adjElem(j),:),edge_nodes(i,:));
            parent(adjElem(j)).locEdge(k1:k2) = repmat(locEdge,1,itrsect(i).num);
            parent(adjElem(j)).coincident(locEdge) = itrsect(i).coincident;
        
            if(isfield(itrsect(i),'nurbsParam'))
                parent(adjElem(j)).nurbsParam(k1:k2) = itrsect(i).nurbsParam;
            end
            if(itrsect(i).multItrsect)
                parent(adjElem(j)).multItrsect = true;
            end
                
        end
    end
end      

for i = 1:size(elem_nodes,1)
    [uniqueNode,ia,ic] = unique(parent(i).nodes);
    parent(i).nItrsect = length(uniqueNode);
    if(parent(i).nItrsect > 1)    
        %  not duplicate nodes
        if (length(parent(i).nodes) == length(uniqueNode))            
            parent(i).nodes = uniqueNode;
            parent(i).edge = parent(i).edge(ia);
            parent(i).locEdge =  parent(i).locEdge(ia);
            parent(i).nSeg = parent(i).nSeg(ia);
            parent(i).seg = parent(i).seg(ia);
            if(isfield(parent(i),'nurbsParam'))
                parent(i).nurbsParam = parent(i).nurbsParam(ia);   
            end
        % remove duplicated nodes which have the same line segments    
        else          
            copy = parent(i);
            parent(i).nodes = uniqueNode;
            parent(i).edge = parent(i).edge(ia);
            parent(i).locEdge = parent(i).locEdge(ia);
            parent(i).nSeg = zeros(1,length(parent(i).nodes));
            parent(i).seg = cell(1,length(parent(i).nodes));
            if (isfield(parent(i),'nurbsParam'))
                parent(i).nurbsParam = cell(1,length(parent(i).nodes));
            end            
            for j = 1:length(copy.nodes)
               parent(i).seg{ic(j)} = [parent(i).seg{ic(j)},copy.seg{j}];
               if(isfield(parent(i),'nurbsParam'))
                    parent(i).nurbsParam{ic(j)} ...
                        = [parent(i).nurbsParam{ic(j)},copy.nurbsParam{j}];
               end
            end
            for j = 1:length(parent(i).nodes)
              [parent(i).seg{j},ia,~] = unique(parent(i).seg{j},'stable');
              if(isfield(parent(i),'nurbsParam'))
                    parent(i).nurbsParam{j} = parent(i).nurbsParam{j}(ia);
              end
              parent(i).nSeg(j) = length(parent(i).seg{j});
            end                    
        end
        
        % detect whether element has 3 intersections per channel, each edge
        % intersect once with channel
        %{
        segs = cat(1,parent(i).seg{:});
        uniqueSegs = unique(segs);
        frequencies = countmember(uniqueSegs,segs);
        if (any(frequencies > 2))
            parent(i).refine = true;           
        end
        %}
    end
end
end