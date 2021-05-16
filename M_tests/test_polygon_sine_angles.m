close all
clear all
path(path,'../../Opt-IGFEM-Curves-2D/M_optimization')
path(path, '../M_geom_toolbox')
rng('shuffle')
%X = [0,1,2,0];
%Y = [0,0,1,1];
%{
X = -1 + 2*rand(4,1);
Y = -1 + 2*rand(4,1);
[X,Y] = poly2ccw(X,Y);
%}
%{
X = [3.272320626598310e-02
     6.983988946080882e-02
     1.287725922939661e-01
     8.826013762007445e-02];
Y = [1.873403754925998e-01
     8.725072238641531e-02
     9.957569102719754e-02
     1.554800467382197e-01];
%}
%{
X = [6.983988946080882e-02
     3.333026915246971e-02
     1.405693903142610e-01
     1.287725922939661e-01];
Y = [8.725072238641531e-02
     6.633203181443051e-03
     3.207070975045748e-02
     9.957569102719754e-02];
%}
%{
channelFile = 'parallel_stephen_start.channel';
readOptions.boundsFile = 'parallel_stephen_start_bounds.channel';
readOptions.generatePolygons = true;
readOptions.nlconfun = @nonlinear_constraints;
readOptions.nlcon.minPolyArea = 1e-4;
readOptions.nlcon.sinMinPolyAngle = sin(0.5*pi/180);
readOptions.mergeTriangles = true;
[channels,~,designParams] = read_channels(channelFile,readOptions);
%}
%polygonFile = 'problemPolygons1';
%load(polygonFile)
%{
figure
plot_channel_polygons(channels.polygons,designParams.vertices2params)
%}

channels.polygons(1).XX = [0,0,1,0;
                           0,0,0,1];

for ii = 1:numel(channels.polygons)
    fprintf('\npolygon number %i \n',ii)
    [sineAngles,dSineAngles] = polygon_sine_angles(channels.polygons(ii).XX(1,:),...
                                                   channels.polygons(ii).XX(2,:));
                                        
    [area,dArea] = polygon_area(channels.polygons(ii).XX(1,:),...
                                channels.polygons(ii).XX(2,:));
    
    nVertices = size(channels.polygons(ii).XX,2);
    nVertices2 = 2*nVertices;
    dSineAnglesFD = nan(nVertices2,nVertices);
    dAreaFD = nan(nVertices2,1);
    del = 1e-7;
    for i = 1:nVertices2
        newXX = channels.polygons(ii).XX;
        newXX(i) = newXX(i) - del;
        newSineAngles1 = polygon_sine_angles(newXX(1,:),newXX(2,:));
        newArea1 = polygon_area(newXX(1,:),newXX(2,:));
        
        newXX = channels.polygons(ii).XX;
        newXX(i) = newXX(i) + del;
        newSineAngles2 = polygon_sine_angles(newXX(1,:),newXX(2,:));
        newArea2 = polygon_area(newXX(1,:),newXX(2,:));
        
        dSineAnglesFD(i,:) = 0.5*(newSineAngles2 - newSineAngles1)/del;
        
        dAreaFD(i) = 0.5*(newArea2 - newArea1 )/del;
        
    end

    dSineAnglesAbsDiff = abs(dSineAngles - dSineAnglesFD);
    dAnglesAbsDiff = abs(dArea - dAreaFD); 
    fprintf('sine angle derivative: max |diff| = %g \n',max(dSineAnglesAbsDiff(:)))

    fprintf('sine angle derivative: max |rel diff| = %g \n',...
             max(dSineAnglesAbsDiff(:)./dSineAngles(:)))
         
    fprintf('area derivative: max |diff| = %g \n',max(dAnglesAbsDiff(:)))     
    fprintf('area derivative: max |rel diff| = %g \n',...
             max(dAnglesAbsDiff(:)./dArea))
end

