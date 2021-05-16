close all
clear
path(path,'../../NURBS/nurbs_toolbox')
geomCstrFile = '../GeomConstrFiles/parallel5x5NE.polygon';
channelFile = '../ChannelFiles/parallel5x5ref.channel';

shift = [0,0];
scaling = 2.5;

[polygons, vertexCoords,...
 vertices2params,...
 restrictedParams.nParams, ...
 restrictedParams.iniVals, ...
 restrictedParams.paramPairs, ...
 sideTriangles] ...
    = read_polygon_file(geomCstrFile );

% scale and shift polygon vertex coordinates
vertexCoords = bsxfun(@plus,scaling*vertexCoords,shift);

% scale and shift restricted parameters initial values
% NOTE: I assume that the restricted parameters are y-coordinates of 
% of control points 
restrictedParams.iniVals = scaling*restrictedParams.iniVals + shift(2);

[channels,designParams] = read_channels(channelFile);

% scale and shift channel end points and NURBS control points
channels.pts = bsxfun(@plus,channels.pts*scaling,shift); 
for i = 1:numel(channels.nurbs)
    channels.nurbs(i).coefs(1:2,:) = bsxfun(@plus,channels.nurbs(i).coefs(1:2,:)*scaling,shift');
end


figure
plot_channel_polygons(polygons, ...
                      vertexCoords, ...
                      channels.nurbs, ...
                      designParams, ...
                      true)  
                  
[filename,pathname] = uiputfile('*.polygon','Save as','../GeomConstrFiles/');
outFile = [pathname,filename];
write_polygon_file(outFile,polygons,vertexCoords,vertices2params, ...
                   restrictedParams.nParams, restrictedParams.iniVals, ...
                   restrictedParams.paramPairs, sideTriangles)
               



