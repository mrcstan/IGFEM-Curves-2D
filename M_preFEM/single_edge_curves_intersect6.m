%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Marcus Tan on 8/20/2013
%%% Last modified date: 12/24/2014
%%% Copyright 2013 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This function finds the intersection of a single edge with multiple
% nurbs curves. 
% INPUT:
% itrsect
% curEdge_nodes
% curEdgeCoords
% currentNode
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
% itrsect.num: number of intersections
% itrsect.nSeg: a length itrsect.num array specifying the number of nurbs
%               curve segments passing through a node 
% itrsect.seg: a itrsect.num cell specifying the nurbs 
%               curve segments passing through a node 
% itrsect.nurbsParam: a itrsect.num cell specfying the location on the
%                     nurbs curve segment that has the intersection point
% itrsect.ipt: a length itrsect.num array indicating the terminal
%                     point of a section of nurbs interface. if the
%                     intersection point is not a terminal point, it is set
%                      to zero.
% itrsect.node: a length itrsect.num array specifying the node number
% itrsect.x
% itrsect.y
% itrsect.potDirichlet 
% if calcItrsectVel == 0
%       itrsect.vx, itrsect.vy: an array of cells of length = itrsect.num.
%           each element of the cell contains a vector of 
%           the velocity of the corresponding node wrt design parameters 
%       itrsect.designParamNum: an array of cells of length = itrsect.num. 
%           each element of the cell contains a vector of
%           the design parameter numbers corresponding to the velocities
%           above
function [itrsect, currentNode,channels] ...
                = single_edge_curves_intersect(itrsect, curEdge_nodes, ...
                                               curEdgeCoords, currentNode, ...
                                               channels, tol, ...
                                               calcItrsectVel)
% loop through all the nurbs curve forming the interface
for j=1:channels.nNurbs % for j       
        % channel is a linear NURBS that can consist of multiple straight
        % segments
        
        
        if (channels.nurbs(j).order == 2)
            [intPts,intPars1,intPars2,intInd,isCoincident] = ...
                        intersect_edges(channels.lineSegs{j}, ...
                                        curEdgeCoords(:)', ...
                                        tol.intersectEdges);
             % convert parametric coordinates from that of a straight segment
             % to that of the linear nurbs
             if (~isempty(intInd))                
                 if(numel(channels.lineSegs{j})>1)
                     inds = channels.nurbs(j).order - 1 + intInd;
                     intPars1 = (1-intPars1).*channels.nurbs(j).knots(inds)' ...
                                +intPars1.*channels.nurbs(j).knots(inds+1)';
                 end
                 intPars1(intPars1 > 1) = 1.0;
                 intPars1(intPars1 < 0) = 0.0;
             end     
            
        else
            [intPts,intPars1,intPars2,isCoincident] = ...
                        curve_edge_nurbs_nurbs_intersect(...
                                channels.nurbs(j).knots,....
                                channels.nurbs(j).coefs(channels.rows{channels.kind(j)},:),...
                                channels.kind(j),...
                                curEdgeCoords,...
                                tol.epsco,tol.epsge);
                                                
        end
         if(any(isCoincident)) % important to maintain itrsect.coincident if it is already true!!!
                itrsect.coincident = true;
         end
       
        for k = 1:length(intPars1)
             % intersection point is one of the end points of
             % the nurbs curve       
             if(abs(intPars1(k))<tol.nurbsParam )                          
                 ipt=channels.contvty(j,1);
             elseif(abs(intPars1(k)-1)<tol.nurbsParam )
                 ipt=channels.contvty(j,2);
             else
                 ipt=0;
             end

             if(ipt>0 && ~ipt_considered(ipt,itrsect))
                 itrsect.potDirichlet = true;
                 [seg1,nSeg1,nurbsParam1]=find_seg_sharing_channels_pt(ipt,channels.contvty);

             elseif(ipt == 0)
                 nSeg1 = 1;
                 seg1 = j;
                 % replace param(1) by a knot value in nurbs.knots if their abs diff is
                 % less than tol
                 %nurbsParam1 = approx_nurbs_param(...
                 %       channels.nurbs(j),closestParam(1),tol.nurbsParam );
                 nurbsParam1 = intPars1(k);
             else
                 break
             end

            % determine if this intersection point is
            % sufficiently far away from the originals nodes
            % so that enrichment node should
            % be added. The nurbs curve is approximated by
            % a straight line passing through the
            % intersection point with a direction equal to
            % the tangent of the curve at that point.
            % This approximately ensures that the adjacent
            % edges also see the same shifted curve
            % when it is shifted to coincide with an
            % original node. Note: this will only work well
            % if halflinewidth is small enough

            %[~,tangent] = nrbdeval(channels.nurbs(j),channels.dnurbs(j),intPars1(k));
            [addEnrichment,locNode,~] = add_enrichment_or_original_node(curEdgeCoords,...
                                                              intPars2(k),...
                                                              intPts(:,k),...
                                                              tol.halfLineWidth);

            % intersection is an enrichment node
            if (addEnrichment)           
                 % check whether the intersection point has been
                 % found on other edges. if so, do not add
                 % enrichment node. if not, add an enrichment node
                 loc = zeros(nSeg1,1);
                 for si = 1:nSeg1
                    [inList,loc(si),channels.nurbsParams{seg1(si)}] ...
                        = find_and_update_sorted_real(nurbsParam1(si),...
                                                      channels.nurbsParams{seg1(si)},...
                                                      tol.nurbsParam);
                 end   
                    
                 if(~inList)                    
                     itrsect.num = itrsect.num+1;
                     num = itrsect.num;
                     itrsect.nSeg(num) = nSeg1;
                     itrsect.seg{num} = seg1;
                     itrsect.ipt(num) = ipt;
                     itrsect.nurbsParam{num} = nurbsParam1;
                     currentNode = currentNode+1;
                     itrsect.node(num) = currentNode;
                     itrsect.x(num) = intPts(1,k);
                     itrsect.y(num) = intPts(2,k); 
                     [~,tangent] = nrbdeval(channels.nurbs(j),channels.dnurbs(j),intPars1(k));
                     itrsect.tangents(:,num) = tangent(1:2);
                     for si = 1:nSeg1
                         ind = seg1(si);
                         if(isempty(channels.enrichNode{ind}))
                             channels.enrichNode{ind} = currentNode;
                         else
                            channels.enrichNode{ind} = [channels.enrichNode{ind}(1:loc(si)-1),...
                                                       currentNode,...
                                                       channels.enrichNode{ind}(loc(si):end)];
                         end
                     end
                 end

             % original node is an intersection on the edge
            else
                 itrsect.num = itrsect.num+1; 
                 num = itrsect.num;
                 itrsect.nSeg(num) = nSeg1;
                 itrsect.seg{num} = seg1;
                 itrsect.ipt(num) = ipt;
                 itrsect.nurbsParam{num} = nurbsParam1;
                 itrsect.x(num) = curEdgeCoords(1,locNode);
                 itrsect.y(num) = curEdgeCoords(2,locNode);
                 itrsect.node(num) = curEdge_nodes(locNode);     
                 [~,tangent] = nrbdeval(channels.nurbs(j),channels.dnurbs(j),intPars1(k));
                 itrsect.tangents(:,num) = tangent(1:2);
                
            end
            if(calcItrsectVel && (~addEnrichment || ~inList))
                 [itrsect.vx{num}, itrsect.vy{num}] ...
                     = intersection_velocities(curEdgeCoords, ...
                                               reshape(channels.lineSegs{j}(intInd(k),:),2,2), ...
                                               tol.vert);
                 itrsect.designParamNum{num} ...
                         = reshape(channels.designParamNum{j}(1:2,(intInd(k):intInd(k)+1)),1,4);                              
            end
        end % k           

end % j
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Marcus Tan on 8/20/2013
%%% Last modified date: 1/5/2014
%%% Copyright 2013 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% this function finds all segments sharing a point
function [seg,nSeg,nurbsParam]=find_seg_sharing_channels_pt(ipt,contvty)
[row,col]=find(contvty==ipt);
nSeg=length(row);
seg=row';
nurbsParam=col'-1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Marcus Tan on 8/20/2013
%%% Last modified date: 8/20/2013
%%% Copyright 2013 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% this function determines whether an intersection point which is a 
% terminal point of a channel/curve has already been found.
function flag=ipt_considered(ipt,itrsect)
flag=false;
if(isfield(itrsect,'ipt'))
    if(any(itrsect.ipt==ipt))
        flag=true;
    end
end
end


%%% Created by Marcus Tan on 1/6/2014
%%% Last modified date: 3/31/2014
%%% Copyright 2014 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% determine if this intersection point is
% sufficiently far away from the originals nodes
% so that enrichment node should be added
function [addEnrichment,locNode,dist] ...
                = add_enrichment_or_original_node(edgeCoords, ...
                                                  edgeItrsectParam, ...
                                                  itrsectPt,tol)

if(edgeItrsectParam < 0.5)
    locNode = 1;
else
    locNode = 2;
end
%line = [itrsectPt',tangent'];
%dist = distancePointLine(edgeCoords(:,locNode)', line); 
dist = norm(edgeCoords(:,locNode)-itrsectPt,2);
if(dist < tol)
    addEnrichment = false;
else
    addEnrichment = true;
    locNode = 0;
end
end


%%% Created by Marcus Tan on 1/5/2014
%%% Last modified date: 1/5/2014
%%% Copyright 2014 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% replace nurbsParam by a knot value in nurbs.knots if their abs diff is
% less than tol
%{
function nurbsParam = approx_nurbs_param(nurbs,nurbsParam,tol)
 nKnot=length(nurbs.knots);
 if(nKnot > 2*nurbs.order)
     for m = nurbs.order+1:nKnot-nurbs.order
         if((nurbs.knots(m)-nurbsParam) > tol)
             break;
         elseif(abs(nurbsParam-nurbs.knots(m)) < tol)
             nurbsParam =nurbs.knots(m);
             break;
         end
     end
 end
end
%}


