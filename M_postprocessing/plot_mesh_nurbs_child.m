%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Marcus Tan on 6/8/2013
%%% Last modified date: 9/3/2012
%%% Copyright 2013 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%
% This function plot original mesh, interface and intersection points between
% mesh and interface and child elements

function []=plot_mesh_nurbs_child(nodeCoord,elemNode,parent,...
                                  showChildNum,showNodeNum,elem2plot)
if (nargin < 4)
    showChildNum = false;
end
if (nargin < 5)
    showNodeNum = false;
end
nElem = size(elemNode,1);
if (nargin < 5 || isempty(elem2plot))
    elem2plot = 1:nElem;
end
ftsize = 20;
ftsize2 = 16;
legendCount = 0;
figure
hold on

hAll = [];
legendStr = [];
nElem2plot = numel(elem2plot);
if (nElem2plot == nElem)
    hElem = trimesh(elemNode,nodeCoord(:,1),nodeCoord(:,2),'color','k');
    hAll=hElem(end);
    legendCount=legendCount+1;
    legendStr{legendCount}='Original mesh';
end

hChild=[];
hLineSource=[];
hInteriorCtrlPt=[];
subd=[30,30];
linespecs.color = 'r';
linespecs.linestyle = '--';
linespecs.width = 2;
%linespecs2={'ko','ks','kd','k^','k<','kv','k>'};
linespecs2={'ro'};
nSpecs=length(linespecs2);
colormap('summer');
cmap=colormap;
rowcmap=size(cmap,1);
mrowcmap=floor(0.75*rowcmap);
for i = 1:nElem2plot
    elemNum = elem2plot(i);
    if (parent(elemNum).type == 2)
        for j = 1:numel(parent(elemNum).child)
            XX = nodeCoord(parent(elemNum).child(j).nodes,:);
            hChild=fill(XX(:,1),XX(:,2),cmap(mrowcmap,:));
              % plot line sources
            for k=1:numel(parent(elemNum).child(j).channelNum)
                XX2 = XX(parent(elemNum).child(j).channelLocNodes(:,k),:);
                hLineSource = plot(XX2(:,1),XX2(:,2),[linespecs.color,linespecs.linestyle],'linewidth',linespecs.width);
            end 
            
            
            
            if (showChildNum)
                XXc = mean(XX,1);
                text(XXc(1),XXc(2),num2str(j),'Color','k',...
                    'fontsize',ftsize2,'BackgroundColor',[.7 .9 .7])
            end
            
            % plot enrichment nodes
            XX2 = XX(parent(elemNum).child(j).locEnNodes,:);
            hInteriorCtrlPt = plot(XX2(:,1),XX2(:,2),linespecs2{1},'MarkerSize',8);
            
            if (showNodeNum)
                for k = 1:numel(parent(elemNum).child(j).nodes)
                    text(XX(k,1),XX(k,2),num2str(parent(elemNum).child(j).nodes(k)),...
                        'Color','k','fontsize',ftsize2)
                end
            end
            
            
        end
    elseif(parent(elemNum).type == 3)
        for j=1:numel(parent(elemNum).child)

            surfspecs.colormap='summer';
            surfspecs.C=j;
            parent(elemNum).child(j).nurbsSurf.coefs(1:2,:,:) = parent(elemNum).child(j).nurbsSurf.coefs(1:2,:,:);
            hChild=nurbs_plot(parent(elemNum).child(j).nurbsSurf,subd,surfspecs);
           
            % plot line sources
            for k=1:numel(parent(elemNum).child(j).channelNum)
                parent(elemNum).child(j).nurbsSeg(k).coefs(1:2,:)=parent(elemNum).child(j).nurbsSeg(k).coefs(1:2,:);
                hLineSource=nurbs_plot(parent(elemNum).child(j).nurbsSeg(k),subd(1),linespecs);
            end       
            % plot interior control points
            nurbsInd=parent(elemNum).child(j).nurbsInd;
            for k=1:length(nurbsInd)
                if(~isempty(nurbsInd{k}))
                    %[sub1,sub2] = ind2sub(parent(elemNum).child(j).nurbsSurf.number,nurbsInd{k});
                    X=parent(elemNum).child(j).nurbsSurf.coefs(1,nurbsInd{k})...
                        ./parent(elemNum).child(j).nurbsSurf.coefs(4,nurbsInd{k});
                    Y=parent(elemNum).child(j).nurbsSurf.coefs(2,nurbsInd{k})...
                        ./parent(elemNum).child(j).nurbsSurf.coefs(4,nurbsInd{k});
                    ind=mod(nurbsInd{k}(1)-1,nSpecs)+1;
                    hInteriorCtrlPt=plot(X,Y,linespecs2{ind},'MarkerSize',8);
                end
            end
            
            %XX = nodeCoord(parent(elemNum).child(j).nodes,:);
            %plot(XX(:,1),XX(:,2),'bd','markersize',8)
            if (showChildNum)
                XX = nodeCoord(parent(elemNum).child(j).nodes,:);
                XX = mean(XX,1);
                text(XX(1),XX(2),num2str(j),'Color','k',...
                    'fontsize',ftsize2,'BackgroundColor',[.7 .9 .7])
            end
            if (showNodeNum)
                XX = nodeCoord(parent(elemNum).child(j).nodes,:);
                for k = 1:numel(parent(elemNum).child(j).nodes)
                    text(XX(k,1),XX(k,2),num2str(parent(elemNum).child(j).nodes(k)),...
                        'Color','k','fontsize',ftsize2)
                end
            end
        end
    end
end
if(~isempty(hChild))
    hAll=[hAll,hChild];
    legendCount=legendCount+1;
    legendStr{legendCount}='Child element';
end
if(~isempty(hLineSource))
    hAll=[hAll,hLineSource];
    legendCount=legendCount+1;
    legendStr{legendCount}='Microchannel';
end

if(~isempty(hInteriorCtrlPt))
    hAll=[hAll,hInteriorCtrlPt];
    legendCount=legendCount+1;
    legendStr{legendCount}='Intersections/Enrichment nodes';
end
xlabel('$\hat{x}$', 'fontsize',ftsize,'Interpreter','latex');
ylabel('$\hat{y}$', 'fontsize',ftsize,'Interpreter','latex');
%
legend(hAll,legendStr)
axis image
set(gca,'fontsize',ftsize)
end

