%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Marcus Tan on 2/6/2015
%%% Copyright 2016 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this function reads a file that describes the polygons used for
% applying geometrical constraints
% FILE FORMAT:
%{
polygons, <number of polygons,p>
v11 v12 ... v1n1
v21 v22 ... v2n2
...
vp1 vp2 ... vpnp
# where v11 v12 etc are vertex numbers

vertexCoords, <number of vertexCoords,c>
x11 x12 
x21 x22
...
xc1 xc2

vertices2params, <number of maps, m>
vp11 vp12
vp21 vp22
...
vpm1 vpm2
# vpij means the design parameter corresponding the j coordinate of the
# ith vertex. If the vertex coordinate does not map to any design
# parameter, the corresponding entry should be nan

# the data from here onwards is only applicable when all polygons are
# triangles
restrictedParams, <number of restricted parameters, r>
val1
val2
...
valr
# restrictedParams are design parameters where the values are
# restricted to be equal to some other design parameters
# these extra parameters are introduced to increase the flexibility of
# the side triangles when applying the geometrical constraints 
# val1, val2, ...are the inital values of these restricted parameters

paramPairs, r
pp11 pp12
pp21 pp22
...
ppr1 ppr2
# (pp11,pp12) is a pair consisting of a restricted parameter and 
# a partner design parameter

sideTriangles, <number of side triangles, s>
triangle1
triangle2
...
trianglep
# these are the side triangles among those listed under polygons
%}
% EXAMPLE FILE (See parallelTwo_NE.polygon):
%{
    polygons,12
    9 8 3 
    10 9 3 
    1 10 2 
    2 10 3 
    7 2 3 
    8 4 3 
    7 3 4 
    12 7 4 
    5 11 4 
    8 5 4 
    11 12 4 
    6 11 5 
    vertexCoords,12
    0 0.19
    0.02 0.19
    0.02 0.1
    0.13 0.1
    0.13 0.01
    0.15 0.01
    0.13 0.19
    0.02 0.01
    0 0.01
    0 0.1
    0.15 0.1
    0.15 0.19
    vertices2params,12
    NaN NaN
    1 2
    5 6
    7 8
    11 12
    NaN NaN
    3 4
    9 10
    NaN 13
    NaN 14
    NaN 15
    NaN 16
    restrictedParams,4
    0.01
    0.1
    0.1
    0.19
    paramPairs,4
    10 13
    6 14
    8 15
    4 16
    sideTriangles,8
    1
    2
    3
    4
    8
    9
    11
    12
%}
function  [polygons,vertexCoords,vertices2params,...
          nRestrictedParams,rParamIniVals,paramPairs,sideTriangles] ...
                = read_polygon_file(inputFile)
polygons = [];
vertexCoords = [];
vertices2params = [];
nRestrictedParams = [];
rParamIniVals = [];
paramPairs = [];
sideTriangles = [];
fileID = fopen(inputFile,'r');
line = fgetl(fileID);
while(ischar(line))
    if (strncmpi(line,'polygons',8))
        splitStr = regexp(line,',','split');
        nPolygons = str2double(splitStr{2});
        polygons = struct('nVertices',0,'connectivity',cell(nPolygons,1));
        for i = 1:nPolygons
            polygons(i).connectivity = str2num(fgetl(fileID));
            polygons(i).nVertices = numel(polygons(i).connectivity);
        end

    elseif (strncmpi(line,'vertexCoords',12))
        splitStr = regexp(line,',','split');
        nVertexCoords = str2double(splitStr{2});
        vertexCoords = fscanf(fileID,'%g %g\n',[2,nVertexCoords])';
    elseif (strncmpi(line,'vertices2params',15))
        splitStr = regexp(line,',','split');
        nVertices2params = str2double(splitStr{2});
        vertices2params = fscanf(fileID,'%g %g\n',[2,nVertices2params])';
    elseif (strncmpi(line,'restrictedParams',16))
        splitStr = regexp(line,',','split');
        nRestrictedParams = str2double(splitStr{2});
        rParamIniVals = fscanf(fileID,'%g \n',[1,nRestrictedParams])';
    elseif (strncmpi(line,'paramPairs',10))
        splitStr = regexp(line,',','split');
        nParamPairs = str2double(splitStr{2});
        paramPairs = fscanf(fileID,'%i %i\n',[2,nParamPairs])';
    elseif (strncmpi(line,'sideTriangles',13))
        splitStr = regexp(line,',','split');
        nSideTriangles = str2double(splitStr{2});
        sideTriangles = fscanf(fileID,'%i\n',[1,nSideTriangles])';
    end
    line = fgetl(fileID);
end
fclose(fileID);
end