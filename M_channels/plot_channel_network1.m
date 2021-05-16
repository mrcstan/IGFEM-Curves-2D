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
ftsize = 20;
% offset nodal values and labels from the nodes so that they have less
% chance of overlapping with other values and labels
offSetNodalLabel = [1,1,0]; % direction
offSetFactor = 0.2;
hold on

nChans = size(chan_nodes,1);
nNodes = max(chan_nodes(:));
lengthScale = inf;
if (isstruct(geom))
    if (length(geom) ~= nChans)
        error('number of nurbs should be equal to the number of rows in chan_nodes')
    end
    specs.width = 2;
    specs.curve = 'r-';
    subd = 100;
    node_coords = zeros(nChans,3);
    chanLength = zeros(nChans,1);
    for i = 1:length(geom)
        nurbs_plot(geom(i),subd,specs);
        endNodes = chan_nodes(i,:);
        node_coords(endNodes,:) = nrbeval(geom(i),[0,1])';
        chanLength(i) = nurbs_arc_length(10,geom(i));
    end
    lengthScale = min(chanLength);

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


showNodes = false;
showChans = false;
showNodalVals = true;
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

if (isstruct(geom))
     for i = 1:nChans
        XX = nrbeval(geom(i),0.5);
        if (showChans)
            text(XX(1),XX(2),XX(3),['<',num2str(i),'>'],'color','r','fontsize',ftsize);
        end
        if (showChanVals)  
            dnurbs = nrbderiv(geom(i));
            [~,tangent] = nrbdeval(geom(i),dnurbs,0.5);
            tangent = 2.0*lengthScale*sign(chanVals(i))*tangent/norm(tangent);
            h = quiver3(XX(1),XX(2),XX(3),tangent(1),tangent(2),tangent(3),0.2,'k-','linewidth',2);
            adjust_quiver_arrowhead_size(h, 10.0);            
            
            text(XX(1),XX(2),XX(3),['(',num2str(i),') ',num2str(chanVals(i),'%15.3g')],...
                'color','r','fontsize',ftsize);
        end    
    end
else
    for i = 1:nChans
        XX = node_coords(chan_nodes(i,:),:); 
        delX = diff(XX);
        normDelX = norm(delX);
        if (normDelX < lengthScale)
            lengthScale = normDelX;
        end
        plot3(XX(:,1),XX(:,2),XX(:,3),'k-')
        if (showChans)                       
            XXm = mean(XX,1);
            text(XXm(1),XXm(2),XXm(3),['<',num2str(i),'>'],'color','r','fontsize',ftsize);
        end
        if (showChanVals)          
            tangent = sign(chanVals(i))*(XX(2,:)-XX(1,:));
            XXm = mean(XX,1);
            h = quiver3(XXm(1),XXm(2),XXm(3),tangent(1),tangent(2),tangent(3),0.2,'k-','linewidth',2);
            adjust_quiver_arrowhead_size(h, 10.0);
            text(XXm(1),XXm(2),XXm(3),['(',num2str(i),') ',num2str(abs(chanVals(i)),'%15.3g')],...
                'color','r','fontsize',ftsize);
        end    
    end
end

offSetNodalLabel = offSetNodalLabel*offSetFactor*lengthScale;
if (showNodes)
    for i = 1:nNodes
        sgn = (-1)^i;
        text(node_coords(i,1)+sgn*offSetNodalLabel(1),...
             node_coords(i,2)+sgn*offSetNodalLabel(2),...
             node_coords(i,3)+sgn*offSetNodalLabel(3),...
             ['<',num2str(i),'>'],'color','b','fontsize',ftsize);
    end
end


if (showNodalVals)
     for i = 1:nNodes
         sgn = (-1)^i;
        text(node_coords(i,1)+sgn*offSetNodalLabel(1),...
             node_coords(i,2)+sgn*offSetNodalLabel(2),...
             node_coords(i,3)+sgn*offSetNodalLabel(3),...
             ['<',num2str(i),'> ',num2str(nodalVals(i),'%15.3g')],...
             'color','b','fontsize',ftsize);
    end
end

xlabel('x', 'fontsize',ftsize);
ylabel('y', 'fontsize',ftsize);
axis image
set(gca,'fontsize',ftsize)

end