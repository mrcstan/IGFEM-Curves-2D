%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Marcus Tan on 7/10/2014
%%% Last modified date: 7/10/2014
%%% Copyright 2014 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this function plots the channel network if the nodal positions 
% or nurbs descriptions of the channels are given
% INPUT:
%   chan_nodes: a number of channels x 2 matrix of the channels and their
%               end nodes
%   geom: an array containing the geom description of each channel. it can be: 
%         i) an array of nurbs structures or
%         ii) number of channels x number of dimensions matrix
%             of the nodal coordinates
%   nodalVals: an array of values at the nodes
%   chanVals: an array of values of the channels
function [] = plot_channel_network(chan_nodes,geom,nodalVals,chanVals)
ftsize = 24;
hold on

nChans = size(chan_nodes,1);
nNodes = max(chan_nodes(:));

if (isstruct(geom))
    if (length(geom) ~= nChans)
        error('number of nurbs should be equal to the number of rows in chan_nodes')
    end
    specs.width = 2;
    specs.curve = 'r-';
    subd = 100;
    node_coords = zeros(nNodes,3);

    for i = 1:length(geom)
        nurbs_plot(geom(i),subd,specs);
        endNodes = chan_nodes(i,:);
        node_coords(endNodes,:) = nrbeval(geom(i),[0,1])';
    end
 

elseif (isnumeric(geom))
    if (size(geom,1) ~= nNodes)
        error('number of rows of geom should equal to the number of nodes')
    end  
    if (size(geom,2) == 2)
        node_coords = [geom,zeros(nNodes,1)];
    elseif(size(geom,2) == 3)
        node_coords = geom;
    else
        error('number of columns of geom should be either 2 or 3')
    end
else
    error('unknown geom variable type')
end

% In order of decreasing precedence, showNodes,showNodalVals
showNodes = false;
showNodalVals = false;
% In order of decreasing precedence, showChans,showChanVals
showChans = false;
showChanVals = true;

if (nargin < 3)
    showNodes = true;
    showChans = true;
end
if (nargin < 3 || isempty(nodalVals))
    showNodalVals = false;
end

if (nargin < 4 || isempty(chanVals))
    showChanVals = false;
end

%scaleChanVal = max(chanVals);
scaleArrow = 20;
if (isstruct(geom))
     for i = 1:nChans
        XX = nrbeval(geom(i),0.5);
        if (showChans)
            text(XX(1),XX(2),XX(3),['(',num2str(i),')'],'color','r','fontsize',ftsize);       
        elseif (showChanVals)  
            dnurbs = nrbderiv(geom(i));
            [~,tangent] = nrbdeval(geom(i),dnurbs,0.5);
            %tangent = 2.0*lengthScale*sign(chanVals(i))*tangent/norm(tangent);
            tangent = sign(chanVals(i))*tangent/norm(tangent)/scaleArrow;
            h = quiver3(XX(1),XX(2),XX(3),tangent(1),tangent(2),tangent(3),0.2,'k-','linewidth',2);
            adjust_quiver_arrowhead_size(h, 10.0);            
            
            text(XX(1),XX(2),XX(3),[num2str(abs(chanVals(i)),'%15.3g')],...
                'color','r','fontsize',ftsize);          
        end
    end
else
    for i = 1:nChans
        XX = node_coords(chan_nodes(i,:),:); 

        plot3(XX(:,1),XX(:,2),XX(:,3),'k-')
        if (showChans)                       
            XXm = mean(XX,1);
            text(XXm(1),XXm(2),XXm(3),['(',num2str(i),')'],'color','r','fontsize',ftsize);
        elseif (showChanVals)          
            tangent = sign(chanVals(i))*(XX(2,:)-XX(1,:));
            tangent = sign(chanVals(i))*tangent/norm(tangent)/scaleArrow;
            XXm = mean(XX,1);
            h = quiver3(XXm(1),XXm(2),XXm(3),tangent(1),tangent(2),tangent(3),0.2,'k-','linewidth',2);
            adjust_quiver_arrowhead_size(h, 10.0);
            text(XXm(1),XXm(2),XXm(3),[num2str(abs(chanVals(i)),'%15.3g')],...
                'color','r','fontsize',ftsize); 
       end    
    end
end

hAlignment = {'center','center'};
nhAlignment = numel(hAlignment);
vAlignment = {'top','bottom'};
nvAlignment = numel(vAlignment);
if (showNodes)
    for i = 1:nNodes
        hAlign = hAlignment{rem((i-1),nhAlignment)+1};
        vAlign = vAlignment{rem((i-1),nvAlignment)+1};
        text(node_coords(i,1),...
             node_coords(i,2),...
             node_coords(i,3),...
             ['<',num2str(i),'>'],...
             'color','b',...
             'fontsize',ftsize,...
             'VerticalAlignment', vAlign,...
             'HorizontalAlignment',hAlign);
    end
elseif (showNodalVals)
     for i = 1:nNodes
        hAlign = hAlignment{rem((i-1),nhAlignment)+1};
        vAlign = vAlignment{rem((i-1),nvAlignment)+1};
        text(node_coords(i,1),...
             node_coords(i,2),...
             node_coords(i,3),...
             ['<',num2str(i),'> ',num2str(nodalVals(i),'%15.3g')],...
             'color','b',...
             'fontsize',ftsize,...
             'VerticalAlignment', vAlign,...
             'HorizontalAlignment',hAlign);
    end
end

set(gcf,'units','normalized','outerposition',[0 0 1 1])
set(gcf,'color',[0.95,0.95,0.95])
set(gca,'fontsize',ftsize)

xlabel('x', 'fontsize',ftsize);
ylabel('y', 'fontsize',ftsize);

axis image 
box on
bounds = [get(gca,'xlim'),get(gca,'ylim')];
dels = [bounds(2)-bounds(1),bounds(4)-bounds(3)];
maxDels = max(dels);
newDels = [max(dels(1),0.1*maxDels),...
           max(dels(2),0.1*maxDels)];

centers = 0.5*[bounds(1)+bounds(2),bounds(3)+bounds(4)];
xlim([centers(1)-0.5*newDels(1),centers(1)+0.5*newDels(1)])
ylim([centers(2)-0.5*newDels(2),centers(2)+0.5*newDels(2)])
%axis off
end