%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Marcus Tan on 8/20/2013
%%% Last modified date: 8/25/2013
%%% Copyright 2013 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%
% This function plot original mesh, interface and intersection points between
% mesh and interface.

function []=plot_mesh_curve_itrsect_junc(node,elem_nodes,channels,itrsect,...
                                         junc,kinks,showElem,showNode,markedElems)
if (nargin < 7)
    showElem = false;
end
if (nargin < 8)
    showNode = false;
end

nElem=size(elem_nodes,1);
if (nargin < 9 || isempty(markedElems))
    colorElems = false;
else
    colorElems = true;
end
ftsize=16;
figure
hold on

legendCount=0;
hJunc = [];
hKinks = [];
hMarkedElem = [];
hElem = trimesh(elem_nodes,node.coords(:,1),node.coords(:,2),'color','k');
hElem = hElem(1);
if (showElem || showNode || colorElems)
    for i=1:nElem 
        XX = node.coords(elem_nodes(i,:),:);
       % plot element numbers
        XXc = mean(XX,1);
        if (colorElems && markedElems(i))
            hMarkedElem = fill(XX(:,1),XX(:,2),'y');
        end
        if(showElem)
            text(XXc(1),XXc(2),num2str(i),'Color','b')
        end

        % node numbers
        %{
        if(showNode)
            text(XX(1,1),XX(1,2),num2str(elem_nodes(i,1)),'Color','k');
            text(XX(2,1),XX(2,2),num2str(elem_nodes(i,2)),'Color','k');
            text(XX(3,1),XX(3,2),num2str(elem_nodes(i,3)),'Color','k');
        end
        %}
    end
end
if (~isempty(junc) || ~isempty(kinks))
    for i = 1:nElem
        if(junc(i).nJunc)
            for j=1:junc(i).nJunc
                %X=channels.pts(junc(i).ipt(j),1);
                %Y=channels.pts(junc(i).ipt(j),2);
                XX = node.coords(junc(i).node(j),:);
                hJunc=plot(XX(:,1),XX(:,2),'bs','markerfacecolor','b','markersize',8);
                if (showNode)
                    text(XX(1),XX(2),num2str(junc(i).node(j)),'color','g')
                end
            end
        end        
        if (kinks(i).nKinks)
            for j = 1:kinks(i).nKinks
                XX = node.coords(kinks(i).node(j),:);
                hKinks=plot(XX(:,1),XX(:,2),'cd','markerfacecolor','c','markersize',8);
                if (showNode)
                    text(XX(1),XX(2),num2str(kinks(i).node(j)),'color','g')
                end
            end
        end
    end
end
%
hAll=hElem;
legendCount=legendCount+1;
legendStr{legendCount}='Mesh';
if(~isempty(hJunc))
    hAll=[hAll,hJunc];
    legendCount=legendCount+1;
    legendStr{legendCount}='Interior branching points';
end

if(~isempty(hKinks))
    hAll=[hAll,hKinks];
    legendCount=legendCount+1;
    legendStr{legendCount}='Interior kinks';
end

if(~isempty(hMarkedElem))
    hAll=[hAll,hMarkedElem];
    legendCount=legendCount+1;
    legendStr{legendCount}='Marked elements';
end

hchannels=[];
linespecs.color = 'r';
linespecs.linestyle = '-';
linespecs.width = 2;
%refWidth = 2;
%refDiam = 7.5e-4;
%magnification = 1;

subd=2000;
for i=1:channels.nNurbs
    %linespecs.width = magnification*refWidth*channels.diams(i)/refDiam; 
    channels.nurbs(i).coefs(1:2,:) = channels.nurbs(i).coefs(1:2,:);
    hchannels=nurbs_plot(channels.nurbs(i),subd,linespecs);
end
if(~isempty(hchannels))
    hAll=[hAll,hchannels];
    legendCount=legendCount+1;
    legendStr{legendCount}='Microchannels';
end

nEdge=length(itrsect);
hItrsect=[];
%
for i=1:nEdge
    if(itrsect(i).num>0)        
        for j=1:itrsect(i).num
            X=itrsect(i).x(j); 
            Y=itrsect(i).y(j);
            hItrsect=plot(X,Y,'ro','markersize',6);
            if(showNode)
                text(X,Y,num2str(itrsect(i).node(j)),'color','r')
            end
        end       
    end  
end
%
%{
for i = node.nOriginalNode+1:node.n_node
    hItrsect = plot(node.coords(i,1),node.coords(i,2),'ro','markersize',8);
    if (showNode)
        text(node.coords(i,1),node.coords(i,2),num2str(i),'color','r')
    end
end
%}
if(~isempty(hItrsect))
    hAll=[hAll,hItrsect];
    legendCount=legendCount+1;
    legendStr{legendCount}='Intersection points';
end

%xlabel('$\hat{x}$', 'fontsize',ftsize,'Interpreter','latex');
%ylabel('$\hat{y}$', 'fontsize',ftsize,'Interpreter','latex');

legend(hAll,legendStr, ...
         3,  'Location', 'NorthOutside');
axis image
set(gca,'fontsize',ftsize)
end

