close all
clear all
path(path,'../../Opt-IGFEM-Curves-2D/M_optimization')
path(path, '../M_geom_toolbox')
X = [0,1,1,0];
Y = [0,0,-1.1,1];

figure
hold on
plot([X,X(1)],[Y,Y(1)],'r-','linewidth',2)
area = polygon_area(X,Y)