%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Marcus Tan on 6/8/2013
%%% Last modified date: 10/10/2013
%%% Copyright 2013 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This function finds the intersection of the edges of the elements with
% the interface consisting of line segments. 
% the equation of a line segment is given by (x,y)=(xx1,yy1)+s(vx,vy)
% where (vx,vy)=(xx2-xx1,yy2-yy1)
% the equation of an edge is given by (x,y)=(x1,y1)+t(wx,wy)
% where (wx,wy)=(x2-x1,y2-y1).
% I assume that when a line segment intersects an element edge, the line
% segment is at least as long as the element edge. The function has not
% been tested when this assumption does not hold. 
% INPUT:
% edge_node
% node.coords
% node.n_node
% itrface.pts
% itrface.contvty
% itrface.junc
% halfLineWidth: half the line segment width within which a terminal node of an edge is
% considered to intersect the line segment
% isConforming: if mesh is conforming stop finding intersections if an edge
% has a maximum of two intersection points.

% OUTPUT:
% itrsect.num: number of intersections
% itrsect.nSeg: a length itrsect.num array specifying the number of line
%               source segments passing through a node 
% itrsect.seg: a itrsect.num x itrsect.nSeg array specifying the line 
%               source segments passing through a node 
% itrsect.itrfacePt: a length itrsect.num array indicating the terminal
%                     point of the line segment. if the intersection point
%                     is not a terminal point, it is set to zero.
% itrsect.node: a length itrsect.num array specifying the node number
% itrsect.x
% itrsect.y
% node.n_node
% node.node_n
% node.coords

function [itrsect,node]=edge_line_segs_intersect(edge_node,node,itrface,halfLineWidth,isConforming)
% stol is the fraction of a line segment below which the intersection point
% is considered different from the line segment terminal points.
%stol=1e-5;
%minso=-stol;
%maxso=1+stol;
% ttol is the fraction of an edge length below which the intersection point
% is considred to be the same as one of the terminal points of the edge.
% increase ttol to reduce sensitivity of defining new enrichment nodes
% ttol=0.001; 
% the angle below which the line segment is considered to be
% parallel to the edge;
% angleTol=1.0*pi/180.0;
% half the line segment width within which a terminal node of an edge is
% considered to intersect the line segment
% halfLineWidth=0.1/40*0.01;
lineWidth=2*halfLineWidth;
nEdge=size(edge_node,1);
nEnrichNode=0;
itrsect=struct('num',cell(nEdge,1));
for j=1:itrface.nSeg
    iPt1=itrface.contvty(j,1);
    iPt2=itrface.contvty(j,2);
    itrface.vx(j)=itrface.pts(iPt2,1)-itrface.pts(iPt1,1);
    itrface.vy(j)=itrface.pts(iPt2,2)-itrface.pts(iPt1,2);
    itrface.length(j)=norm([itrface.vx(j),itrface.vy(j)],2);
    % stol is the fraction of a line segment below which the intersection point
    % is considered different from the line segment terminal points.
    itrface.stol(j)=halfLineWidth/itrface.length(j);
    itrface.minso(j)=-itrface.stol(j);
    itrface.maxso(j)=1+itrface.stol(j);
end
for i=1:nEdge % for i
    INd1=edge_node(i,1);
    INd2=edge_node(i,2);
    x1=node.coords(INd1,1);
    x2=node.coords(INd2,1);
    y1=node.coords(INd1,2);
    y2=node.coords(INd2,2);
    itrsect(i).num=0;
    itrsect(i).nSeg=[];
    itrsect(i).seg=[];
    itrsect(i).pt=[];
    itrsect(i).node=[];
    itrsect(i).x=[];
    itrsect(i).y=[];
    edgeLen=norm([x2-x1,y2-y1],2);
    for j=1:itrface.nSeg % for j
        iPt1=itrface.contvty(j,1);
        xx1=itrface.pts(iPt1,1);
        yy1=itrface.pts(iPt1,2);
        iPt2=itrface.contvty(j,2);
        xx2=itrface.pts(iPt2,1);
        yy2=itrface.pts(iPt2,2);
        xxL=min(xx1,xx2)-halfLineWidth;
        xxR=max(xx1,xx2)+halfLineWidth;
        yyL=min(yy1,yy2)-halfLineWidth;
        yyR=max(yy1,yy2)+halfLineWidth;
        if(~((x1<xxL && x2<xxL) || (y1<yyL && y2<yyL) || ...
            (x1>xxR && x2>xxR)  || (y1>yyR && y2>yyR))) % if1
            vx=itrface.vx(j);
            vy=itrface.vy(j);
            segLen=itrface.length(j);
            wx=x2-x1;
            wy=y2-y1;       
            d1=xx1-x1;
            d2=yy1-y1;
            A=[-vx,wx;-vy,wy];
            % angle between the direction vector of the line segment
            % and the edge
            cosAngle=dot([vx,vy],[wx,wy])/(segLen*edgeLen);
            % the angle below which the line segment is considered to be
            % parallel to the edge;
            cosAngleTol=sqrt(1-(halfLineWidth/edgeLen)^2);
            stol=itrface.stol(j);
            minso=itrface.minso(j);
            maxso=itrface.maxso(j);
            % line segment not parallel to the element edge
            if(abs(cosAngle)<cosAngleTol) % end if det  
                % solve system of equation formed by equating the line
                % segment and the edge equations.
                soln=A\[d1;d2];
                so=soln(1);
                to=soln(2);              
                % line segment intersects the edge
                if(so>=minso && so<=maxso)
                     % calculate perpendicular vectors and distances from edge...
                    % nodes to line segments
                    [sNode2seg,dist_sNode2seg,eNode2seg,dist_eNode2seg]=...
                      node2seg(edge_node(i,:),node.coords,j,itrface);
                     % line segment intersects the edge
                    if(dot(sNode2seg,eNode2seg)<0)
                        nSeg0=1;
                        seg0=j;
                        pt0=0;
                        % check if the intersection is at a terminal point of a
                        % line segment
                        if(abs(so)<stol)
                            startPt=itrface.contvty(j,1);
                            startPtFlag=pt_considered(startPt,itrsect(i));
                            if(startPtFlag)
                                continue
                            else
                                [seg0,nSeg0,pt0]=find_seg_sharing_pt(startPt,...
                                                  itrface.contvty);             
                            end
                        elseif(abs(so-1)<stol)
                            endPt=itrface.contvty(j,2);
                            endPtFlag=pt_considered(endPt,itrsect(i));
                            if(endPtFlag)
                                continue
                            else
                                [seg0,nSeg0,pt0]=find_seg_sharing_pt(endPt,...
                                                  itrface.contvty);   
                            end
                        end          
                    
                        % if the distance of an edge node is within a
                        % perpendicular distance of halfLineWidth from the line
                        % segment, then the intersection is taken to be the
                        % edge node
                        if(dist_sNode2seg<halfLineWidth)
                            sNode=edge_node(i,1);
                            itrsect(i).num=itrsect(i).num+1;
                            num=itrsect(i).num;
                            itrsect(i).nSeg(num)=nSeg0;
                            itrsect(i).seg{num}=seg0;
                            itrsect(i).pt(num)=pt0;
                            itrsect(i).node(num)=sNode;
                            itrsect(i).x(num)=node.coords(sNode,1);
                            itrsect(i).y(num)=node.coords(sNode,2);                       
                        elseif(dist_eNode2seg<halfLineWidth)
                            eNode=edge_node(i,2);
                            itrsect(i).num=itrsect(i).num+1;
                            num=itrsect(i).num;
                            itrsect(i).nSeg(num)=nSeg0;
                            itrsect(i).seg{num}=seg0;
                            itrsect(i).pt(num)=pt0;
                            itrsect(i).node(num)=eNode;
                            itrsect(i).x(num)=node.coords(eNode,1);
                            itrsect(i).y(num)=node.coords(eNode,2);   

                        else
                            nEnrichNode=nEnrichNode+1;
                            currentNode=node.n_node+nEnrichNode;
                            coord=[xx1+so*vx,yy1+so*vy];
                            itrsect(i)=update_itrsect(itrsect(i),...
                                seg0,nSeg0,pt0,currentNode,coord);  
                            node.node_n(currentNode)=currentNode;
                            node.coords(currentNode,1:2)=coord;
                        end
                    elseif(dist_sNode2seg<halfLineWidth)
                        nSeg0=1;
                        seg0=j;
                        pt0=0;
                        % check if the intersection is at a terminal point of a
                        % line segment
                        if(abs(so)<stol)
                            startPt=itrface.contvty(j,1);
                            startPtFlag=pt_considered(startPt,itrsect(i));
                            if(startPtFlag)
                                continue
                            else
                                [seg0,nSeg0,pt0]=find_seg_sharing_pt(startPt,...
                                                  itrface.contvty);             
                            end
                        elseif(abs(so-1)<stol)
                            endPt=itrface.contvty(j,2);
                            endPtFlag=pt_considered(endPt,itrsect(i));
                            if(endPtFlag)
                                continue
                            else
                                [seg0,nSeg0,pt0]=find_seg_sharing_pt(endPt,...
                                                  itrface.contvty);   
                            end
                        end        
                        sNode=edge_node(i,1);
                        itrsect(i).num=itrsect(i).num+1;
                        num=itrsect(i).num;
                        itrsect(i).nSeg(num)=nSeg0;
                        itrsect(i).seg{num}=seg0;
                        itrsect(i).pt(num)=pt0;
                        itrsect(i).node(num)=sNode;
                        itrsect(i).x(num)=node.coords(sNode,1);
                        itrsect(i).y(num)=node.coords(sNode,2); 
                    elseif(dist_eNode2seg<halfLineWidth)
                        nSeg0=1;
                        seg0=j;
                        pt0=0;
                        % check if the intersection is at a terminal point of a
                        % line segment
                        if(abs(so)<stol)
                            startPt=itrface.contvty(j,1);
                            startPtFlag=pt_considered(startPt,itrsect(i));
                            if(startPtFlag)
                                continue
                            else
                                [seg0,nSeg0,pt0]=find_seg_sharing_pt(startPt,...
                                                  itrface.contvty);             
                            end
                        elseif(abs(so-1)<stol)
                            endPt=itrface.contvty(j,2);
                            endPtFlag=pt_considered(endPt,itrsect(i));
                            if(endPtFlag)
                                continue
                            else
                                [seg0,nSeg0,pt0]=find_seg_sharing_pt(endPt,...
                                                  itrface.contvty);   
                            end
                        end  
                        eNode=edge_node(i,2);
                        itrsect(i).num=itrsect(i).num+1;
                        num=itrsect(i).num;
                        itrsect(i).nSeg(num)=nSeg0;
                        itrsect(i).seg{num}=seg0;
                        itrsect(i).pt(num)=pt0;
                        itrsect(i).node(num)=eNode;
                        itrsect(i).x(num)=node.coords(eNode,1);
                        itrsect(i).y(num)=node.coords(eNode,2); 
                    else
                        continue
                    end
                end
            % line segment parallel to the element edge    
            else
                sNode=edge_node(i,1);
                eNode=edge_node(i,2);
                % calculate perpendicular vectors and distances from edge...
                % nodes to line segments
                [sNode2seg,dist_sNode2seg,eNode2seg,dist_eNode2seg]=...
                      node2seg(edge_node(i,:),node.coords,j,itrface);
                % edge is considered to be aligned with the line segment if
                % it lies within a rectangular region of width
                % 2xhalfLineWidth
                if(dot(sNode2seg,eNode2seg)<0 || ...
                   (dist_sNode2seg<halfLineWidth && ...
                    dist_eNode2seg<halfLineWidth))
                    % project edge nodes onto the line segment 
                    sNodeCoord=node.coords(sNode,:)+sNode2seg;
                    eNodeCoord=node.coords(eNode,:)+eNode2seg;
                    wwx=eNodeCoord(1)-sNodeCoord(1);
                    wwy=eNodeCoord(2)-sNodeCoord(2);
                    dd1=xx1-sNodeCoord(1);
                    dd2=yy1-sNodeCoord(2);
                    if (abs(vx)>lineWidth)
                        so=-dd1/vx;
                        s1=(wwx-dd1)/vx;
                    elseif (abs(vy)>lineWidth)
                        so=-dd2/vy;
                        s1=(wwy-dd2)/vy;
                    else
                        continue
                    end
                else
                    continue
                end
                startPt=itrface.contvty(j,1);
                endPt=itrface.contvty(j,2);
                startPtFlag=pt_considered(startPt,itrsect(i));                            
                endPtFlag=pt_considered(endPt,itrsect(i));
                % t=0 node of edge intersects with line segment
                if(so>=minso && so<=maxso) % if so
                    % both nodes of the edge lie on the line segment
                    if(s1>=minso && s1<=maxso) % if s1a
                        if(abs(so)<stol && abs(s1-1)<stol)
                            if(startPtFlag && endPtFlag)
                                continue;
                            elseif(startPtFlag)
                                [seg1,nSeg1,pt1]=find_seg_sharing_pt(endPt,...
                                              itrface.contvty);
                                 itrsect(i)=update_itrsect(itrsect(i),...
                                            seg1,nSeg1,pt1,eNode,...
                                            node.coords(eNode,:)); 
                            elseif(endPtFlag)
                                [seg0,nSeg0,pt0]=find_seg_sharing_pt(startPt,...
                                              itrface.contvty);
                                itrsect(i)=update_itrsect(itrsect(i),...
                                            seg0,nSeg0,pt0,sNode,...
                                            node.coords(sNode,:));
                            end
                        elseif(abs(so-1)<stol && abs(s1)<stol)
                            if(startPtFlag && endPtFlag)
                                continue;
                            elseif(startPtFlag)
                                [seg0,nSeg0,pt0]=find_seg_sharing_pt(endPt,...
                                              itrface.contvty);
                                itrsect(i)=update_itrsect(itrsect(i),...
                                            seg0,nSeg0,pt0,sNode,...
                                            node.coords(sNode,:)); 
                            elseif(endPtFlag)
                                [seg1,nSeg1,pt1]=find_seg_sharing_pt(startPt,...
                                              itrface.contvty);
                                 itrsect(i)=update_itrsect(itrsect(i),...
                                    seg1,nSeg1,pt1,eNode,...
                                    node.coords(eNode,:));           
                            end
                        elseif(abs(so)<stol)
                            if(startPtFlag)
                                itrsect(i)=update_itrsect(itrsect(i),...
                                            j,1,0,eNode,...
                                            node.coords(eNode,:));
                            else
                                [seg0,nSeg0,pt0]=find_seg_sharing_pt(startPt,...
                                              itrface.contvty);
                                itrsect(i)=update_itrsect(itrsect(i),seg0,nSeg0,pt0,...
                                        sNode,node.coords(sNode,:));
                                itrsect(i)=update_itrsect(itrsect(i),j,1,0,...
                                        eNode,node.coords(eNode,:));
                            end
                        elseif(abs(s1)<stol)
                            if(startPtFlag)
                               itrsect(i)=update_itrsect(itrsect(i),...
                                            j,1,0,sNode,...
                                            node.coords(sNode,:));
                            else
                               [seg1,nSeg1,pt1]=find_seg_sharing_pt(startPt,...
                                              itrface.contvty);
                               itrsect(i)=update_itrsect(itrsect(i),j,1,0,...
                                        sNode,node.coords(sNode,:));
                               itrsect(i)=update_itrsect(itrsect(i),seg1,nSeg1,pt1,...
                                        eNode,node.coords(eNode,:));
                            end
                        elseif(abs(so-1)<stol)
                            if(endPtFlag)
                                itrsect(i)=update_itrsect(itrsect(i),...
                                            j,1,0,eNode,...
                                            node.coords(eNode,:));
                            else
                               [seg0,nSeg0,pt0]=find_seg_sharing_pt(endPt,...
                                              itrface.contvty);
                               itrsect(i)=update_itrsect(itrsect(i),seg0,nSeg0,pt0,...
                                        sNode,node.coords(sNode,:));
                               itrsect(i)=update_itrsect(itrsect(i),j,1,0,...
                                        eNode,node.coords(eNode,:));
                            end
                        elseif(abs(s1-1)<stol)
                            if(endPtFlag)
                                 itrsect(i)=update_itrsect(itrsect(i),...
                                            j,1,0,sNode,...
                                            node.coords(sNode,:));
                            else
                               [seg1,nSeg1,pt1]=find_seg_sharing_pt(endPt,...
                                              itrface.contvty);
                               itrsect(i)=update_itrsect(itrsect(i),j,1,0,...
                                        sNode,node.coords(sNode,:));
                               itrsect(i)=update_itrsect(itrsect(i),seg1,nSeg1,pt1,...
                                        eNode,node.coords(eNode,:));
                            end
                        else                      
                            itrsect(i)=update_itrsect(itrsect(i),j,1,0,...
                                        sNode,node.coords(sNode,:));
                            itrsect(i)=update_itrsect(itrsect(i),j,1,0,...
                                        eNode,node.coords(eNode,:));
                        end                                          
                    % an end point A of the line segment terminates 
                    % on the element edge
                    else
                        if(abs(so)<stol)
                            if(~startPtFlag)
                                [seg0,nSeg0,pt0]=find_seg_sharing_pt(startPt,...
                                              itrface.contvty);
                                itrsect(i)=update_itrsect(itrsect(i),...
                                            seg0,nSeg0,pt0,sNode,node.coords(sNode,:));
                            end                                
                        elseif(abs(so-1)<stol)
                            if(~endPtFlag)
                                [seg0,nSeg0,pt0]=find_seg_sharing_pt(endPt,...
                                              itrface.contvty);
                                itrsect(i)=update_itrsect(itrsect(i),...
                                            seg0,nSeg0,pt0,sNode,node.coords(sNode,:));
                            end
                                
                        end
                        % A is the s=1 end of the
                        % line segment if s1>maxso 
                        % A is different from the
                        % t=1 node if abs(so-1)>tol
                        if(s1>maxso && abs(so-1)>stol) % if s1b
                            if(endPtFlag)
                                 continue;
                            else
                                 [seg1,nSeg1,pt1]=find_seg_sharing_pt(...
                                     endPt,itrface.contvty);
                            end
                            nEnrichNode=nEnrichNode+1;
                            currentNode=node.n_node+nEnrichNode;
                            coord=[xx1+vx,yy1+vy];
                            itrsect(i)=update_itrsect(itrsect(i),...
                                            seg1,nSeg1,pt1,currentNode,coord);  
                            node.node_n(currentNode)=currentNode;
                            node.coords(currentNode,1:2)=coord;
                        % A is the s=0 end of the
                        % line segment    
                        % A is different from the
                        % t=0 node if abs(so)>tol
                        elseif(s1<minso && abs(so)>stol)
                            if(startPtFlag)
                                continue;
                            else
                                 itrsect(i).num=itrsect(i).num+1;
                                 num=itrsect(i).num;
                                 [seg1,nSeg1,pt1]=find_seg_sharing_pt(...
                                    itrface.contvty(j,1),...
                                    itrface.contvty);
                            end
                            nEnrichNode=nEnrichNode+1;
                            currentNode=node.n_node+nEnrichNode;
                            coord=[xx1,yy1];
                            itrsect(i)=update_itrsect(itrsect(i),...
                                            seg1,nSeg1,pt1,currentNode,coord);   
                            node.node_n(currentNode)=currentNode;
                            node.coords(currentNode,1:2)=coord;
                        end % if s1b
                    end % if s1a
                % t=1 node of edge intersects with line segment and 
                % an end point B of the line segment terminates 
                % on the edge of the element
                elseif(s1>=minso && s1<=maxso)
                        if(abs(s1)<stol)
                            if(~startPtFlag)
                                [seg1,nSeg1,pt1]=find_seg_sharing_pt(startPt,...
                                              itrface.contvty);
                                itrsect(i)=update_itrsect(itrsect(i),...
                                            seg1,nSeg1,pt1,eNode,node.coords(eNode,:));
                            end                                
                        elseif(abs(s1-1)<stol)
                            if(~endPtFlag)
                                [seg1,nSeg1,pt1]=find_seg_sharing_pt(endPt,...
                                              itrface.contvty);
                                itrsect(i)=update_itrsect(itrsect(i),...
                                            seg1,nSeg1,pt1,eNode,node.coords(eNode,:));
                            end
                                
                        end
                        % B is the s=1 end of the
                        % line segment if so>maxso
                        % B is different from the
                        % t=1 node if abs(s1-1)>tol
                        if(so>maxso && abs(s1-1)>stol)
                            if(endPtFlag)
                                continue;
                            else
                                 itrsect(i).num=itrsect(i).num+1;
                                 num=itrsect(i).num;
                                 [seg1,nSeg1,pt1]=find_seg_sharing_pt(...
                                     endPt,itrface.contvty);
                            end
                            nEnrichNode=nEnrichNode+1;
                            currentNode=node.n_node+nEnrichNode;
                            coord=[xx1+vx,yy1+vy];
                            itrsect(i)=update_itrsect(itrsect(i),...
                                            seg1,nSeg1,pt1,currentNode,coord);
                            node.node_n(currentNode)=currentNode;
                            node.coords(currentNode,1:2)=coord;
                        % B is the s=0 end of the
                        % line segment if so<minso  
                        % B is different from the
                        % t=0 node if abs(s1)>tol
                        elseif(so<minso && abs(s1)>stol)
                            if(startPtFlag)
                                continue;                
                            else
                                itrsect(i).num=itrsect(i).num+1;
                                num=itrsect(i).num;
                                [seg1,nSeg1,pt1]=find_seg_sharing_pt(...
                                         startPt,itrface.contvty);
                            end
                            nEnrichNode=nEnrichNode+1;
                            currentNode=node.n_node+nEnrichNode;
                            coord=[xx1,yy1];
                            itrsect(i)=update_itrsect(itrsect(i),...
                                           seg1,nSeg1,pt1,currentNode,coord);
                            
                            node.node_n(currentNode)=currentNode;
                            node.coords(currentNode,1:2)=coord;
                        end  
                end % if so           
            end % end if det
        end % end if1 
        for k = 1:itrsect(i).num
            itrsect(i).seg{k} = unique(itrsect(i).seg{k});
            itrsect(i).nSeg(k) = length(itrsect(i).seg{k});
        end
        if(isConforming && itrsect(i).num>2)
            break
        end
    end
end

node.n_node=size(node.coords,1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Marcus Tan on 7/23/2013
%%% Last modified date: 7/23/2013
%%% Copyright 2013 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% this function finds all segments sharing a point
function [seg,nSeg,pt]=find_seg_sharing_pt(ipt,contvty)
[row,col]=find(contvty==ipt);
nSeg=length(row);
seg=row';
if(nSeg>1)
    pt=ipt;
else
    pt=0;
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Marcus Tan on 7/24/2013
%%% Last modified date: 7/24/2013
%%% Copyright 2013 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% this function determines whether an intersection point which is a 
% terminal point of a line segment has already been found.
function flag=pt_considered(pt,itrsect)
flag=false;
if(isfield(itrsect,'pt'))
    if(any(itrsect.pt==pt))
        flag=true;
    end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Marcus Tan on 7/25/2013
%%% Last modified date: 7/25/2013
%%% Copyright 2013 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function itrsectNew=update_itrsect(itrsectOld,seg,nSeg,pt,node,nodeCoords)
    itrsectNew=itrsectOld;
    itrsectNew.num=itrsectNew.num+1;
    num=itrsectNew.num;
    itrsectNew.nSeg(num)=nSeg;
    itrsectNew.seg{num}=seg;
    itrsectNew.pt(num)=pt;
    itrsectNew.node(num)=node;
    itrsectNew.x(num)=nodeCoords(1);
    itrsectNew.y(num)=nodeCoords(2);   
end




