close all
clear all
path(path,'../M_channels')
path(path,'../channelFiles')
path(path,'../../NURBS/nurbs_toolbox')
path(path,'../../Opt-IGFEM-Curves-2D/M_optimization')
path(path, '../M_geom_toolbox')

channelFile = 'parallel5x5ref.channel';
readOptions.boundsFile = [];
readOptions.sampleFile = [];
readOptions.sampleNum = 14;
readOptions.polygonFile = 'parallel5x5NE.polygon';
readOptions.polygonFig = [];
readOptions.nlconfun = @nonlinear_constraints;
readOptions.nlcon.minPolyArea = 0.001*0.15*0.2;
readOptions.nlcon.sinMinPolyAngle = sin(0.5*pi/180);
readOptions.mergeTriangles = false;
[channels,~,designParams,restrictedParams] = preprocess_channels(channelFile,readOptions);
ylim([0,0.2])
set(gcf,'units','normalized','outerposition',[0 0 1 1])
set(gcf,'color','white')
set(gca,'units','normalized','position',[0.12,0.15,0.75,0.65])
axis off
box on
rectangle('Position',[0,0,0.15,0.2],'edgeColor','k','clipping','off')

figure
plot_channel_polygons(channels.polygons, ...
                      channels.vertexCoords, ...
                      channels.nurbs, ...
                      designParams, ...
                      true)
set(gcf,'units','normalized','outerposition',[0 0 1 1])
set(gcf,'color','white')
set(gca,'units','normalized','position',[0.12,0.15,0.75,0.65])
axis off
box on
rectangle('Position',[0,0,0.15,0.2],'edgeColor','k','clipping','off')

[Aeq,beq] = linear_equality_constraints(restrictedParams,...
                                        designParams.nParams);                  
[lb,ub] = params_bounds(designParams,restrictedParams);