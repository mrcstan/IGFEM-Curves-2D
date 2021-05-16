close all
clear all
path(path,'../ChannelFiles')
path(path, '../M_channels')
path(path, '../M_geom_toolbox')
path(path, '../../NURBS/nurbs_toolbox')
path(path,'../../vflib/mex')
path(path,'../../MatlabUsefulFunctions/export_fig')
path(path,'../../MatlabUsefulFunctions/GrTheory')
ftsz = 30;

channelFile = 'parallel2_start.channel'; 


[channels,designParams] = read_channels(channelFile);
nChannels = size(channels.contvty,1);

[vertex2channel,smallEdgeList,nOrgVertices,nVertices] ...
                = channel2vertex_vertex2channel_matrix(channels,'v2c');

% find original vertex coordinates
channel2vertex = channel2vertex_vertex2channel_matrix(channels,'c2v');
orgVertexCoord = nan(nVertices,2);
orgVertexCoord(1:nOrgVertices,:) = channels.pts;
for i= nOrgVertices+1:nVertices
    chanNum = find(channel2vertex(:,i));
    ctrlPtNum = channel2vertex(chanNum,i);
    orgVertexCoord(i,:) = channels.nurbs(chanNum).coefs(1:2,ctrlPtNum)';
end
%{
desiredNedges = 10;
smallEdgeList = [1,2;
             2,3;
             2,4;
             3,5;
             4,5;
             4,6;
             5,7;
             6,7;
             7,8];
%}
rng('shuffle')
vertexCoord = [0,0.19;[0.14*rand(nVertices-2,1)+0.005,0.19*rand(nVertices-2,1)+0.005];0.15,0.01];
vertexCoord([nOrgVertices,nVertices],:) = vertexCoord([nVertices,nOrgVertices],:); 

constraint = [];

figure
set(gcf,'units','normalized','outerposition',[0 0 1 1])
set(gcf,'color',[0.95,0.95,0.95])
hold on
plot(vertexCoord(:,1),vertexCoord(:,2),'bd','markersize',12,'markerfacecolor','b')


axis image
axis off
rectangle('Position',[0,0,0.15,0.2],'edgeColor','k','clipping','off')
%set(gca,'units','normalized','position',[0.05,0.05,0.8,0.8])
set(gca,'fontsize',30)

DT = delaunayTriangulation(vertexCoord);

E = edges(DT);



figure
set(gcf,'units','normalized','outerposition',[0 0 1 1])
set(gcf,'color',[0.95,0.95,0.95])
hold on
trimesh(DT.ConnectivityList,vertexCoord(:,1),vertexCoord(:,2),'color','k','linestyle','-');
plot(vertexCoord(:,1),vertexCoord(:,2),'bd','markersize',12,'markerfacecolor','b')


[small2largeCol1,small2largeCol2,matchNumber] ...
    = mx_subgraph_monomorphism(smallEdgeList',nVertices,E',nVertices,[1,nOrgVertices]);
if (matchNumber > 0)
    fprintf('match satisfying constraint found \n');
elseif (matchNumber < 0)
    fprintf('match NOT satisfying constraint found \n');
else
    fprintf('NO match found \n');
end

axis image
axis off
rectangle('Position',[0,0,0.15,0.2],'edgeColor','k','clipping','off')
%set(gca,'units','normalized','position',[0.05,0.05,0.8,0.8])
set(gca,'fontsize',30)

if (matchNumber)
    figure
    set(gcf,'units','normalized','outerposition',[0 0 1 1])
    set(gcf,'color',[0.95,0.95,0.95])
    hold on
    trimesh(DT.ConnectivityList,vertexCoord(:,1),vertexCoord(:,2),'color','k','linestyle','-');

    vertexCoord(small2largeCol1,:) = vertexCoord(small2largeCol2,:); 
    for i = 1:size(smallEdgeList,1)
        edgeXX =  vertexCoord(smallEdgeList(i,:),:);
        plot(edgeXX(:,1),edgeXX(:,2),'r-','linewidth',3)
    end
    plot(vertexCoord(:,1),vertexCoord(:,2),'bd','markersize',12,'markerfacecolor','b')
    
    axis image
    axis off
    rectangle('Position',[0,0,0.15,0.2],'edgeColor','k','clipping','off')
    %set(gca,'units','normalized','position',[0.05,0.05,0.8,0.8])
    set(gca,'fontsize',30)
    % check if the orientation (ccw/cw) of the vertices are the same as before
    cycle = grCycleBasis(smallEdgeList);
    cycle1Edges = smallEdgeList(logical(cycle(:,1)),:); 
    G = sparse(cycle1Edges(:,1),cycle1Edges(:,2),1,nVertices,nVertices);
    discVertices = graphtraverse(G,cycle1Edges(1,1));
    
    if (ispolycw(orgVertexCoord(discVertices,1),orgVertexCoord(discVertices,2)) ...
            ~= ispolycw(vertexCoord(discVertices,1),vertexCoord(discVertices,2)))
        warning('vertices have opposite orientation compared to original')
    else
        figure
        set(gcf,'units','normalized','outerposition',[0 0 1 1])
        set(gcf,'color',[0.95,0.95,0.95])
        hold on
        for i = 1:size(vertexCoord,1)
            plot(vertexCoord(i,1),vertexCoord(i,2),'bd','markersize',12,'markerfacecolor','b')
        end
        for i = 1:size(smallEdgeList,1)
            edgeXX =  vertexCoord(smallEdgeList(i,:),:);
            plot(edgeXX(:,1),edgeXX(:,2),'r-','linewidth',3)
        end
        axis image
        axis off
        rectangle('Position',[0,0,0.15,0.2],'edgeColor','k','clipping','off')
        %set(gca,'units','normalized','position',[0.05,0.05,0.8,0.8])
        set(gca,'fontsize',30)
    end
end
