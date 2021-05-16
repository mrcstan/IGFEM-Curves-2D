%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Marcus Tan on 8/20/2013
%%% Last modified date: 5/8/2014
%%% Copyright 2014 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%
% This function plot original mesh, interface and intersection points between
% mesh and interface.

function []=plot_mesh_labels(node_coords,elem_nodes,showElem,showNode,edge_nodes,edgeMarkers)
ftsize=20;
figure
hold on
nElem=size(elem_nodes,1);
%

trimesh(elem_nodes,node_coords(:,1),node_coords(:,2),'color','k')
if (nargin< 3)
	return
end
if (nargin< 4)
	showNode = false;
end
if (nargin< 5)
    edge_nodes = [];
    edgeMarkers = [];
end

if (nargin < 6)
    edgeMarkers = false(size(edge_nodes,1),1);
end

if (showElem || showNode)
    for i=1:nElem 
        XX = node_coords(elem_nodes(i,:),:);
       % plot element numbers
        XXc = mean(XX,1);
        XXmark = 0.9*XX(1,:)+0.05*XX(2,:)+0.05*XX(3,:); 
        if(showElem)
            text(XXc(1),XXc(2),num2str(i),'Color','b','fontsize',ftsize)
            plot(XXmark(1),XXmark(2),'db')
        end
        % node numbers
        if(showNode)
            text(XX(1,1),XX(1,2),num2str(elem_nodes(i,1)),'Color','k','fontsize',ftsize);
            text(XX(2,1),XX(2,2),num2str(elem_nodes(i,2)),'Color','k','fontsize',ftsize);
            text(XX(3,1),XX(3,2),num2str(elem_nodes(i,3)),'Color','k','fontsize',ftsize);
        end
        %
    end
end
for i = 1:size(edge_nodes,1)
    XXc = mean(node_coords(edge_nodes(i,:),:),1); 
    if(edgeMarkers(i))
        text(XXc(1),XXc(2),num2str(i),'Color','r','fontsize',ftsize);
    else
        text(XXc(1),XXc(2),num2str(i),'Color','g','fontsize',ftsize);
    end
end

%

xlabel('x/L', 'fontsize',ftsize);
ylabel('y/L', 'fontsize',ftsize);

%legend(hAll,legendStr, ...
%         3,  'Location', 'NorthOutside');
axis image
set(gca,'fontsize',ftsize)

end

