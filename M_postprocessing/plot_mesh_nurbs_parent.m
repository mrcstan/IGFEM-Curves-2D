%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Marcus Tan on 6/8/2013
%%% Last modified date: 7/7/2014
%%% Copyright 2013 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%
% This function plot original mesh, interface and intersection points between
% mesh and interface and parent elements

function []=plot_mesh_nurbs_parent(nodeCoord,...
                                   nOriginalNode,...
                                   elemNode,...
                                   parent,...
                                   channels,...
                                   showElemNum,...
                                   showNodeNum,...
                                   elem2plot)
if (nargin < 6)
    showElemNum = false;
end
if (nargin < 7)
    showNodeNum = false;
end

ftsize = 20;
ftsize2 = 16;
legendCount=0;
figure
hold on
nElem=size(elemNode,1);
if (nargin < 8 || isempty(elem2plot))
    elem2plot = 1:nElem;
end

hAll = [];
legendStr = [];
nElem2plot = numel(elem2plot);
if (nElem2plot == nElem)
    hElem = trimesh(elemNode,nodeCoord(:,1),nodeCoord(:,2),'color','k');
    hAll=hElem(end);
    legendCount=legendCount+1;
    legendStr{legendCount}='Original mesh';
end

hParent=[];
hHasChannel=[];
hItrsectEnrichNode=[];
hItrsectNode=[];
hnurbsseg=[];
hSurf=[];
colormap('bone');
cmap=colormap;
rowcmap=size(cmap,1);
mrowcmap=floor(0.75*rowcmap);
linespecs.color = 'b';
linespecs.linestyle = '--';
linespecs.width=2;
subd=50; % number of pts to plot for each nurb curve segment
for i = 1:nElem2plot
    elemNum = elem2plot(i);
    if(parent(elemNum).type > 1)
        XX = nodeCoord(parent(elemNum).nodes,:);
        hParent=fill(XX(1:3,1),XX(1:3,2),cmap(mrowcmap,:));
        % parent element number
        if(showElemNum)
            XXc = mean(XX(1:3,:),1);
            text(XXc(1),XXc(2),num2str(elemNum),'Color','k',...
                'fontsize',ftsize2,'BackgroundColor',[.7 .9 .7]);
        end
        % node numbers
        if(showNodeNum)
            for j = 1:numel(parent(elemNum).nodes)
               text(XX(j,1),XX(j,2),num2str(parent(elemNum).nodes(j)),...
                   'Color','k','fontsize',ftsize2);
           end
        end
        % plot enrichment nodes and intersection nodes
        for j=1:numel(parent(elemNum).nodes)
            if(parent(elemNum).nodes(j) > nOriginalNode)
                hItrsectEnrichNode=plot(XX(j,1),XX(j,2),'ko','linewidth',1,'MarkerSize',8);
            elseif(parent(elemNum).nSeg(j) > 0)
                hItrsectNode=plot(XX(j,1),XX(j,2),'r^','MarkerSize',6);
            end   
        end 
        
        if (parent(elemNum).type == 2)
            for j = 1:numel(parent(elemNum).channelNum)
                hnurbsseg = plot(XX(parent(elemNum).channelLocNodes(:,j),1),...
                                 XX(parent(elemNum).channelLocNodes(:,j),2),...
                                 [linespecs.color,linespecs.linestyle],'linewidth',linespecs.width);
            end
        end
        % plot nurbs curve segments
        if (parent(i).type == 3)
            for j=1:numel(parent(elemNum).channelNum)
                hnurbsseg = nurbs_curve_n_tangent_plot(parent(elemNum).nurbsSeg(j),subd,linespecs);
                % hLineSource=nurbs_plot(parent(i).nurbsSeg(j),subd,linespecs);
            end
        end
        %
    elseif(parent(i).type == 1)
        XX = nodeCoord(parent(elemNum).nodes,:);
        hHasChannel = trimesh([1,2,3],XX(1:3,1),XX(1:3,2),'color','c');
        hHasChannel =  hHasChannel(1);
         % element number
        if(showElemNum)
            XXc = mean(XX(1:3,:),1);
            text(XXc(1),XXc(2),num2str(elemNum),'Color','k',...
                'fontsize',ftsize2,'BackgroundColor',[.7 .9 .7]);
        end
        % node numbers
        if(showNodeNum)
            for j = 1:numel(parent(elemNum).nodes)
               text(XX(j,1),XX(j,2),num2str(parent(elemNum).nodes(j)),...
                   'fontsize',ftsize2,'Color','k');
           end
        end
        % plot intersection nodes
        %
        for j=1:numel(parent(elemNum).nodes)
            if(parent(elemNum).nSeg(j)>0)
                hItrsectNode=plot(XX(j,1),XX(j,2),'r^','MarkerSize',6);
            end
        end
        %
        % plot line sources
        for j=1:numel(parent(elemNum).channelNum)
            X = nodeCoord(parent(elemNum).channelNodes(:,j),1);
            Y = nodeCoord(parent(elemNum).channelNodes(:,j),2);
            hnurbsseg =plot(X,Y,'b','linewidth',2);
            h = quiver(X(1),Y(1),X(2)-X(1),Y(2)-Y(1),0.5,'k--','linewidth',3);
            adjust_quiver_arrowhead_size(h, 10.0);
        end
    end
end
hchannels=[];

%
subd=1000;
linespecs.color = 'r';
linespecs.linestyle = '-';
linespecs.width=2;

for i=1:channels.nNurbs
    hchannels=nurbs_plot(channels.nurbs(i),subd,linespecs);
end
%

if(~isempty(hParent))
    hAll=[hAll,hParent];
    legendCount=legendCount+1;
    legendStr{legendCount}= 'Parent element';
end
if(~isempty(hHasChannel))
    hAll=[hAll, hHasChannel];
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
if(~isempty(hnurbsseg))
    hAll=[hAll,hnurbsseg];
      legendCount=legendCount+1;
    legendStr{legendCount}='NURBS segments';
end
if(~isempty(hSurf))
    hAll=[hAll,hSurf];
      legendCount=legendCount+1;
    legendStr{legendCount}='Parent element NURBS surface';
end
if(~isempty(hchannels))
    hAll=[hAll,hchannels];
    legendCount=legendCount+1;
    legendStr{legendCount}='Original NURBS Interface';
end
xlabel('x', 'fontsize',ftsize,'fontweight','b');
ylabel('y', 'fontsize',ftsize,'fontweight','b');

legend(hAll,legendStr)
axis image
set(gca,'fontsize',ftsize)
end

