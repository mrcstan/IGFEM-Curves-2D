%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Marcus Tan on 6/8/2013
%%% Last modified date: 9/3/2012
%%% Copyright 2013 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%
% This function plot original mesh, interface and intersection points between
% mesh and interface and child elements

function []=plot_mesh_nurbs_child(nodeCoord,elemNode,parent,showChild)
if (nargin < 4)
    showChild = false;
end
ftsize=16;
legendCount=0;
figure
hold on
nElem=size(elemNode,1);

hElem = trimesh(elemNode,nodeCoord(:,1),nodeCoord(:,2),'color','k');
hAll=hElem(end);
legendCount=legendCount+1;
legendStr{legendCount}='Background mesh';


hChild=[];
hLineSource=[];
hInteriorCtrlPt=[];
subd=[30,30];
linespecs.curve='r--';
linespecs.width=2;
%linespecs2={'ko','ks','kd','k^','k<','kv','k>'};
linespecs2={'ro'};
nSpecs=length(linespecs2);
for i=1:nElem
    if(parent(i).isParent)
        for j=1:parent(i).nChild
            %
    
            surfspecs.colormap='summer';
            surfspecs.C=j;
            parent(i).child(j).nurbsSurf.coefs(1:2,:,:) = parent(i).child(j).nurbsSurf.coefs(1:2,:,:);
            hChild=nurbs_plot(parent(i).child(j).nurbsSurf,subd,surfspecs);
            %
            % plot line sources
            for k=1:parent(i).child(j).nLineSource
                parent(i).child(j).nurbsSeg(k).coefs(1:2,:)=parent(i).child(j).nurbsSeg(k).coefs(1:2,:);
                hLineSource=nurbs_plot(parent(i).child(j).nurbsSeg(k),subd(1),linespecs);
            end       
            % plot interior control points
            %
            nurbsInd=parent(i).child(j).nurbsInd;
            for k=1:length(nurbsInd)
                if(~isempty(nurbsInd{k}))
                    X=parent(i).child(j).nurbsSurf.coefs(1,nurbsInd{k}(1),nurbsInd{k}(2))...
                        /parent(i).child(j).nurbsSurf.coefs(4,nurbsInd{k}(1),nurbsInd{k}(2));
                    Y=parent(i).child(j).nurbsSurf.coefs(2,nurbsInd{k}(1),nurbsInd{k}(2))...
                        /parent(i).child(j).nurbsSurf.coefs(4,nurbsInd{k}(1),nurbsInd{k}(2));
                    ind=mod(nurbsInd{k}(1)-1,nSpecs)+1;
                    hInteriorCtrlPt=plot(X,Y,linespecs2{ind},'MarkerSize',8);
                end
            end
            %
            XX = nodeCoord(parent(i).child(j).node,:);
            plot(XX(:,1),XX(:,2),'bd','markersize',8)
            if (showChild)
                XX = nodeCoord(parent(i).child(j).node,:);
                XX = mean(XX,1);
                text(XX(1),XX(2),num2str(j),'Color','k','BackgroundColor',[.7 .9 .7])
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

