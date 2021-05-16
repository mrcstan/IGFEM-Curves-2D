path(path, '../M_geom_toolbox')
path(path, '../M_FEM')
close all
xv = [0,rand,0]';
yv = [0,0,rand]';
xmin = min(xv);
xmax = max(xv);
ymin = min(yv);
ymax = max(yv);
%xpts = 0.5;
%ypts = 9e-6;
xpts = xmin+(xmax-xmin)*rand(3,1);
ypts = ymin+(ymax-ymin)*rand(3,1);


dd = 1e-5;
trimesh([1,2,3],xv,yv)
hold on
for i = 1:numel(xpts)
    plot(xpts(i),ypts(i),'ro','markersize',14,'markerfacecolor','r')
    text(xpts(i),ypts(i),num2str(i),'fontsize',40)
end
xlabel('x')
ylabel('y')

[in,locFaces,Xloc] = inTriangle(xpts,ypts,xv,yv,dd)
