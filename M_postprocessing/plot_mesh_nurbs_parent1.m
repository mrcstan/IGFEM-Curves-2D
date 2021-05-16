%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Marcus Tan on 6/8/2013
%%% Last modified date: 7/7/2014
%%% Copyright 2013 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%
% This function plot original mesh, interface and intersection points between
% mesh and interface and parent elements

function []=plot_mesh_nurbs_parent(nodeCoord,nOriginalNode,elemNode,parent,itrface,showElem,showNode)
if (nargin < 6)
    showElem = false;
end
if (nargin < 7)
    showNode = false;
end

ftsize=16;
legendCount=0;
figure
hold on
nElem=size(elemNode,1);

hElem = trimesh(elemNode,nodeCoord(:,1),nodeCoord(:,2),'color','k');
hAll=hElem(end);
legendCount=legendCount+1;
legendStr{legendCount}='Original mesh';


hParent=[];
hHasLineSource=[];
hItrsectEnrichNode=[];
hItrsectNode=[];
hLineSource=[];
hSurf=[];
colormap('bone');
cmap=colormap;
rowcmap=size(cmap,1);
mrowcmap=floor(0.75*rowcmap);
linespecs.curve='b--';
linespecs.width=2;
subd=50; % number of pts to plot for each nurb curve segment
for i=1:nElem
    if(parent(i).isParent)
        XX = nodeCoord(parent(i).node,:);
        hParent=fill(XX(1:3,1),XX(1:3,2),cmap(mrowcmap,:));
        % parent element number
        if(showElem)
            XXc = mean(XX(1:3,:),1);
            text(XXc(1),XXc(2),num2str(i),'Color','k','BackgroundColor',[.7 .9 .7]);
        end
        % node numbers
        if(showNode)
            for j = 1:parent(i).nNode
               text(XX(j,1),XX(j,2),num2str(parent(i).node(j)),'Color','k');
           end
        end
        % plot enrichment nodes and intersection nodes
        for j=1:parent(i).nNode
            if(parent(i).node(j) > nOriginalNode)
                hItrsectEnrichNode=plot(XX(j,1),XX(j,2),'ko','linewidth',1,'MarkerSize',8);
            elseif(parent(i).nSeg(j)>0)
                hItrsectNode=plot(XX(j,1),XX(j,2),'r^','MarkerSize',6);
            end   
        end 
        % plot nurbs curve segments
        %
        for j=1:parent(i).nLineSource
            hLineSource=nurbs_curve_n_tangent_plot(parent(i).nurbsSeg(j),subd,linespecs);
            % hLineSource=nurbs_plot(parent(i).nurbsSeg(j),subd,linespecs);
        end
        %
    elseif(parent(i).hasLineSource)
        XX = nodeCoord(parent(i).node,:);
        hHasLineSource = trimesh([1,2,3],XX(1:3,1),XX(1:3,2),'color','c');
        hHasLineSource = hHasLineSource(1);
         % element number
        if(showElem)
            XXc = mean(XX(1:3,:),1);
            text(XXc(1),XXc(2),num2str(i),'Color','k','BackgroundColor',[.7 .9 .7]);
        end
        % node numbers
        if(showNode)
            for j = 1:parent(i).nNode
               text(XX(j,1),XX(j,2),num2str(parent(i).node(j)),'Color','k');
           end
        end
        % plot intersection nodes
        %
        for j=1:parent(i).nNode
            if(parent(i).nSeg(j)>0)
                hItrsectNode=plot(XX(j,1),XX(j,2),'r^','MarkerSize',6);
            end
        end
        %
        % plot line sources
        for j=1:parent(i).nLineSource
            X = nodeCoord(parent(i).lineSource(j,:),1);
            Y = nodeCoord(parent(i).lineSource(j,:),2);
            hLineSource=plot(X,Y,'b','linewidth',2);
            h = quiver(X(1),Y(1),X(2)-X(1),Y(2)-Y(1),0.5,'k--','linewidth',3);
            adjust_quiver_arrowhead_size(h, 10.0);
        end
    end
end
hItrface=[];
%
subd=1000;
linespecs.curve='r';
linespecs.width=2;

for i=1:itrface.nNurbs
    hItrface=nurbs_plot(itrface.nurbs(i),subd,linespecs);
end
%
%
if(~isempty(hParent))
    hAll=[hAll,hParent];
    legendCount=legendCount+1;
    legendStr{legendCount}= 'Parent element';
end
if(~isempty(hHasLineSource))
    hAll=[hAll,hHasLineSource];
    legendCount=legendCount+1;
    legendStr{legendCount}= 'Element with line source';
end
if(~isempty(hItrsectEnrichNode))
    hAll=[hAll,hItrsectEnrichNode];
    legendCount=legendCount+1;
    legendStr{legendCount}='Enrichment nodes';
end
if(~isempty(hItrsectNode))
    hAll=[hAll,hItrsectNode];
      legendCount=legendCount+1;
    legendStr{legendCount}='Intersection nodes';
end
if(~isempty(hLineSource))
    hAll=[hAll,hLineSource];
      legendCount=legendCount+1;
    legendStr{legendCount}='NURBS segments';
end
if(~isempty(hSurf))
    hAll=[hAll,hSurf];
      legendCount=legendCount+1;
    legendStr{legendCount}='Parent element NURBS surface';
end
if(~isempty(hItrface))
    hAll=[hAll,hItrface];
    legendCount=legendCount+1;
    legendStr{legendCount}='Original NURBS Interface';
end
xlabel('x', 'fontsize',ftsize,'fontweight','b');
ylabel('y', 'fontsize',ftsize,'fontweight','b');

legend(hAll,legendStr)
axis image
set(gca,'fontsize',ftsize)
end

