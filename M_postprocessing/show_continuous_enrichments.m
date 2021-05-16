%% Using matlab quad2d function for integration of nurbs regions
% Gauss - Quadrature for NURBS element

% Preliminary settings
close all
% clear

format long e;

% Adding required toolbox pathes and reading input file
path(path, '../NURBS/nurbs_toolbox') % to add nurbs toolbox path
path(path, '../MatlabUsefulFunctions/export_fig')
ftsz = 34;
subd=[30,30];
subdCurve = 100;
specs.colormap = 'jet';
lnspecs.curve = 'w-';
lnspecs.width = 2;
%
elem_num = 12;
nurbsSurf1=mesh.elem.parent(elem_num).child(1).nurbsSurf;
nurbsSurf2=mesh.elem.parent(elem_num).child(2).nurbsSurf;
nurbsSeg = mesh.elem.parent(elem_num).nurbsSeg;
iPEN=4:mesh.elem.parent(elem_num).nNode;


specs.C = 0;
specs.surf = false;
figure(1)
nodes = mesh.elem.elem_node(elem_num,:);
trimesh([1,2,3],mesh.node.coords(nodes,1),mesh.node.coords(nodes,2),'linestyle','-','color','k','linewidth',2);
hold on
nurbsSeg.coefs(3,2) = 1.0;
nurbs_plot(nurbsSeg, subdCurve, lnspecs)

figure(2)
hold on
trimesh([1,2,3],mesh.node.coords(nodes,1),mesh.node.coords(nodes,2),'linestyle','-','color','k','linewidth',2);
nurbsSeg.coefs(3,2) = 0.0;
nurbsSeg.coefs(3,1) = 1.0;
nurbs_plot(nurbsSeg, subdCurve, lnspecs);

figure(3)
hold on 
trimesh([1,2,3],mesh.node.coords(nodes,1),mesh.node.coords(nodes,2),'linestyle','-','color','k','linewidth',2);
nurbsSeg.coefs(3,1) = 0.0;
lnspecs.curve = 'b-';
nurbs_plot(nurbsSeg, subdCurve, lnspecs);

specs.surf = true;

[~,iCENP,iCEN]=intersect(mesh.elem.parent(elem_num).node(iPEN),...
                        mesh.elem.parent(elem_num).child(1).node,'stable');
nurbsInd=mesh.elem.parent(elem_num).child(1).nurbsInd(iCEN);
[nurbsShape1,~]=create_nurbs_enrichment_function(nurbsSurf1,nurbsInd);
figure(1)
nurbs_plot(nurbsShape1(3),subd, specs);

figure(2)
nurbs_plot(nurbsShape1(2),subd, specs);


[~,iCENP,iCEN]=intersect(mesh.elem.parent(elem_num).node(iPEN),...
                        mesh.elem.parent(elem_num).child(2).node,'stable');
nurbsInd=mesh.elem.parent(elem_num).child(2).nurbsInd(iCEN);
nShape = length(nurbsInd);
[nurbsShape2,~]=create_nurbs_enrichment_function(nurbsSurf2,nurbsInd);
figure(1)
nurbs_plot(nurbsShape2(3),subd, specs);

figure(2)
nurbs_plot(nurbsShape2(2),subd, specs);

figure(1)
grid off
axis off
%caxis([0, 1])
h = colorbar;
set(h,'fontsize',ftsz)
view(35,55)
%

elem_num = 15;
nurbsSurf3=mesh.elem.parent(elem_num).child(1).nurbsSurf;
nurbsSurf4=mesh.elem.parent(elem_num).child(2).nurbsSurf;
nurbsSeg2 = mesh.elem.parent(elem_num).nurbsSeg;
iPEN = 4:mesh.elem.parent(elem_num).nNode;

figure(2)
nodes = mesh.elem.elem_node(elem_num,:);
trimesh([1,2,3],mesh.node.coords(nodes,1),mesh.node.coords(nodes,2),'linestyle','-','color','k','linewidth',2);
nurbsSeg2.coefs(3,3) = 1.0;
lnspecs.curve = 'w-';
nurbs_plot(nurbsSeg2, subdCurve, lnspecs);

figure(1)
nodes = mesh.elem.elem_node(elem_num,:);
trimesh([1,2,3],mesh.node.coords(nodes,1),mesh.node.coords(nodes,2),'linestyle','-','color','k','linewidth',2);
nurbsSeg2.coefs(3,3) = 0.0;
lnspecs.curve = 'w-';
specs.C = 0;
nurbs_plot(nurbsSeg2, subdCurve, lnspecs);
nurbs_plot(nurbsSurf3, subd, specs);
nurbs_plot(nurbsSurf4, subd, specs);

[~,iCENP,iCEN]=intersect(mesh.elem.parent(elem_num).node(iPEN),...
                        mesh.elem.parent(elem_num).child(1).node,'stable');
nurbsInd = mesh.elem.parent(elem_num).child(1).nurbsInd(iCEN);
%nShape = length(nurbsInd);
[nurbsShape3,~]=create_nurbs_enrichment_function(nurbsSurf3,nurbsInd);

figure(2)
nurbs_plot(nurbsShape3(1),subd, specs);

[~,iCENP,iCEN]=intersect(mesh.elem.parent(elem_num).node(iPEN),...
                        mesh.elem.parent(elem_num).child(2).node,'stable');
nurbsInd = mesh.elem.parent(elem_num).child(2).nurbsInd(iCEN);
%nShape = length(nurbsInd);
[nurbsShape4,~]=create_nurbs_enrichment_function(nurbsSurf4,nurbsInd);

nurbs_plot(nurbsShape4(1),subd, specs);

grid off
axis off
caxis([0, 1])
h = colorbar;
set(h,'fontsize',ftsz)
view(40,55)

figure(3)
trimesh([1,2,3],mesh.node.coords(nodes,1),mesh.node.coords(nodes,2),'linestyle','-','color','k','linewidth',2);
nurbsSeg2.coefs(3,3) = 0.0;
lnspecs.curve = 'b-';
nurbs_plot(nurbsSeg2, subdCurve, lnspecs);
view(35,55)
grid off
axis off