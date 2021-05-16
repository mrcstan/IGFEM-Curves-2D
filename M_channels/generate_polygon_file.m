% See read_polygon_file.m for file format and example
close all
clear all
path(path,'../channelFiles')
path(path,'../../NURBS/nurbs_toolbox')
path(path, '../M_geom_toolbox')
path(path, '../M_preFEM')
path(path,'../../MatlabUsefulFunctions/export_fig')
channelFile = 'parallel2x2ref.channel';

options.mergeTriangles = false;
options.orientation = 'NE';

[channels,designParams] = read_channels(channelFile);
[designParams,channels.designParamNum] ...
        = design_params2channel_params_map(designParams,...
                                           channels);

[polygons,vertexCoords,vertices2params,...
          nRestrictedParams,rParamIniVals,paramPairs,sideTriangles] ...
                = channel_polygons(channels,designParams.nParams,options);

figure
plot_channel_polygons(polygons, ...
                      vertexCoords, ...
                      channels.nurbs, ...
                      designParams, ...
                      true)   

[filename,pathname] = uiputfile('*.polygon','Save as','../GeomConstrFiles/');
outFile = [pathname,filename];

write_polygon_file(outFile, polygons,vertexCoords,vertices2params, ...
                   nRestrictedParams, rParamIniVals, paramPairs, sideTriangles)
