%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Marcus Tan on 8/20/2013
%%% Last modified date: 5/8/2014
%%% Copyright 2014 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%
% This function plot original mesh, interface and intersection points between
% mesh and interface.

function plot_mesh_curve(nodeCoords,elemNodes,elemVals,channels,options)
ftsize=30;
hold on

legendCount=0;

if nargin < 5
    options = struct;
end
if ~isfield(options,'showMesh')
    options.showMesh = true;
end
if ~isfield(options,'showElemNums')
    options.showElemNums = false;
end

if ~isfield(options,'showElemVals')
    options.showElemVals = true;
end

if ~isfield(options,'showNodeNums')
    options.showNodeNums = false;
end

if ~isfield(options,'showDiamVals')
    options.showDiamVals = false;
end

if ~isfield(options,'showCtrlPts')
    options.showCtrlPts = true;
end

if ~isfield(options,'bounds')
    options.bounds = [min(nodeCoords(:,1)),max(nodeCoords(:,1)),...
                      min(nodeCoords(:,2)),max(nodeCoords(:,2))];
end

hAll = [];

if ~isempty(elemVals) && options.showElemVals
    if numel(elemVals) ~= size(elemNodes,1)
        error('number of element values must equal to number of elements')
    end
    p = patch('Faces',elemNodes,'Vertices',nodeCoords,'FaceColor','b');
    clear cdata
    minVal = min(elemVals);
    maxVal = max(elemVals);
    if minVal ~= maxVal
        set(gca,'CLim',[minVal,maxVal])
    end
    set(p, 'FaceColor','flat','CData',elemVals,'CDataMapping','scaled')
end

if options.showMesh
    hElem = [];
    nElem=size(elemNodes,1);
    for i=1:nElem
        X=[nodeCoords(elemNodes(i,1),1); nodeCoords(elemNodes(i,2),1); ...
           nodeCoords(elemNodes(i,3),1); ...
           nodeCoords(elemNodes(i,1),1)];
        Y=[nodeCoords(elemNodes(i,1),2); nodeCoords(elemNodes(i,2),2); ...
           nodeCoords(elemNodes(i,3),2); ...
           nodeCoords(elemNodes(i,1),2)];
        hElem=plot(X,Y,'color',[0.7,0.7,0.7],'LineWidth',1);
        if(options.showElemNums)
            text(mean(X),mean(Y),num2str(i),'Color','k',...
                'BackgroundColor',[.7 .9 .7],'fontsize',ftsize)
        end
        % node numbers
        if(options.showNodeNums)
            text(X(1),Y(1),num2str(elemNodes(i,1)),'Color','k','fontsize',ftsize);
            text(X(2),Y(2),num2str(elemNodes(i,2)),'Color','k','fontsize',ftsize);
            text(X(3),Y(3),num2str(elemNodes(i,3)),'Color','k','fontsize',ftsize);
        end
        %
    end
    hAll= [hAll,hElem];
    legendCount = legendCount+1;
    legendStr{legendCount}='Mesh';
end




linespecs.linestyle = '-';
linespecs.width = 3;
%minDiam = 3.5e-4;
%maxDiam = 1.5e-3;
minDiam = 0.0;
maxDiam = 5.6e-4;
delDiam = maxDiam - minDiam;
colormap jet;
cm = colormap; % returns the current color map
nColors = size(cm,1);

subd=2000;
colorIDs = min(nColors,max(1,floor((channels.diams-minDiam)/delDiam * nColors)));
if (nnz(unique(colorIDs)) == 1)
    useColorID = false;
    linespecs.color = 'r';
else
    useColorID = true;
end

hChannel = [];
hCtrlpts = [];

for i = 1:channels.nNurbs
    %linespecs.width = magnification*refWidth*channels.diams(i)/refDiam; 
    if (useColorID)
        linespecs.color = cm(colorIDs(i),:);
    end
    hChannel=nurbs_plot(channels.nurbs(i),subd,linespecs);
    
    if (options.showDiamVals)
        XX = nrbeval(channels.nurbs(i),0.5);
        % diameter in mm
        text(XX(1),XX(2),XX(3),[num2str(channels.diams(i)*1000,3)], ...
            'color','k','fontsize',ftsize);
    end
    if (options.showCtrlPts)
        for j = 1:channels.nurbs(i).number
            hCtrlpts = plot(channels.nurbs(i).coefs(1,j)/channels.nurbs(i).coefs(4,j),...
                            channels.nurbs(i).coefs(2,j)/channels.nurbs(i).coefs(4,j),...
                            'bd','markerfacecolor','b','markersize',12);
            %(channels.nurbs(i).coefs(1,j),channels.nurbs(i).coefs(2,j),...
            %     num2str(j-1),'fontsize',ftsize)
        end
    end
end
%{
caxis([minDiam,maxDiam])
colorbar
%}

%{
if(~isempty(hChannel))
    hAll=[hAll,hChannel];
    legendCount=legendCount+1;
    legendStr{legendCount}='Microchannel';
end

if(~isempty(hCtrlpts))
    hAll=[hAll,hCtrlpts];
    legendCount=legendCount+1;
    legendStr{legendCount}='Control points';
end


legend(hAll,legendStr, ...
         3,  'Location', 'NorthOutside');
%}

set(gcf,'color',[0.95,0.95,0.95])
set(gcf,'units','normalized','outerposition',[0 0 1 1])
set(gca,'units','normalized','position',[0.15,0.15,0.7,0.7])

rectangle('Position', ...
        [options.bounds(1),options.bounds(3),...
         options.bounds(2)-options.bounds(1),...
         options.bounds(4)-options.bounds(3)],'linewidth',1,...
            'edgeColor','k','clipping','off'); 
        
axis image

set(gca,'fontsize',ftsize)
xlim(options.bounds(1:2))
ylim(options.bounds(3:4))

axis off

end

