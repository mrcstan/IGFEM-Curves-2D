function write_polygon_file(outFile, ...
                            polygons,...
                            vertexCoords,...
                            vertices2params,...
                            nRestrictedParams,...
                            rParamIniVals,...
                            paramPairs,...
                            sideTriangles)
    
fileID = fopen(outFile,'w');

fprintf(fileID,'polygons,%i\n',numel(polygons));
for i = 1:numel(polygons)
    fprintf(fileID,'%i ',polygons(i).connectivity);
    fprintf(fileID,'\n');
end
fprintf(fileID,'vertexCoords,%i\n',size(vertexCoords,1));
for i = 1:size(vertexCoords,1)
    fprintf(fileID,'%g %g\n',vertexCoords(i,:));
end
fprintf(fileID,'vertices2params,%i\n',size(vertices2params,1));
for i = 1:size(vertices2params,1)
    fprintf(fileID,'%i %i\n',vertices2params(i,:));
end
fprintf(fileID,'restrictedParams,%i\n',nRestrictedParams);
fprintf(fileID,'%g\n',rParamIniVals);
fprintf(fileID,'paramPairs,%i\n',size(paramPairs,1));
for i = 1:size(paramPairs,1)
    fprintf(fileID,'%i %i\n',paramPairs(i,:));
end
fprintf(fileID,'sideTriangles,%i\n',numel(sideTriangles));
fprintf(fileID,'%i\n',sideTriangles);
fclose(fileID);
end