%% Using matlab quad2d function for integration of nurbs regions
% Gauss - Quadrature for NURBS element

% Preliminary settings
% close all
% clear

format long e;

% Adding required toolbox pathes and reading input file
path(path, '../NURBS/nurbs_toolbox') % to add nurbs toolbox path

subd=[30,30];
specs.colormap='default';

%
elem_num=1;
nurbsSurf1=mesh.elem.parent(elem_num).child(1).nurbsSurf;
nurbsSurf2=mesh.elem.parent(elem_num).child(2).nurbsSurf;
nurbsSurf3=mesh.elem.parent(elem_num).child(3).nurbsSurf;
nurbsSurf4=mesh.elem.parent(elem_num).child(4).nurbsSurf;
iPEN=4:mesh.elem.parent(elem_num).nNode;

%{
specs.C=1;
nurbs_plot(nurbsSurf, subd);
specs.C=2;
nurbs_plot(nurbsSurf2, subd);
specs.C=3;
nurbs_plot(nurbsSurf3, subd);
specs.C=4;
nurbs_plot(nurbsSurf4, subd,specs);
%}
%

figure
hold on
[~,iCENP,iCEN]=intersect(mesh.elem.parent(elem_num).node(iPEN),...
                        mesh.elem.parent(elem_num).child(1).node,'stable');
nurbsInd=mesh.elem.parent(elem_num).child(1).nurbsInd(iCEN);
nShape=length(nurbsInd)
[nurbsShape1,~]=create_nurbs_enrichment_function(nurbsSurf,nurbsInd);
for i=1:2
    %figure
    nurbs_plot(nurbsSurf1, subd);
    %hold on
    nurbs_plot(nurbsShape1(i),subd)
end



[~,iCENP,iCEN]=intersect(mesh.elem.parent(elem_num).node(iPEN),...
                        mesh.elem.parent(elem_num).child(2).node,'stable');
nurbsInd=mesh.elem.parent(elem_num).child(2).nurbsInd(iCEN);
nShape=length(nurbsInd)
[nurbsShape2,~]=create_nurbs_enrichment_function(nurbsSurf2,nurbsInd);
for i=1:2
    %figure
    nurbs_plot(nurbsSurf2, subd);
    %hold on
    nurbs_plot(nurbsShape2(i),subd)
end

[~,iCENP,iCEN]=intersect(mesh.elem.parent(elem_num).node(iPEN),...
                        mesh.elem.parent(elem_num).child(3).node,'stable');
nurbsInd=mesh.elem.parent(elem_num).child(3).nurbsInd(iCEN);
nShape=length(nurbsInd)
[nurbsShape3,~]=create_nurbs_enrichment_function(nurbsSurf3,nurbsInd);
for i=1:3
    %figure
    nurbs_plot(nurbsSurf3, subd);
    %hold on
    nurbs_plot(nurbsShape3(i),subd)
end

[~,iCENP,iCEN]=intersect(mesh.elem.parent(elem_num).node(iPEN),...
                        mesh.elem.parent(elem_num).child(4).node,'stable');
nurbsInd=mesh.elem.parent(elem_num).child(4).nurbsInd(iCEN);
nShape=length(nurbsInd)
[nurbsShape4,~]=create_nurbs_enrichment_function(nurbsSurf4,nurbsInd);
for i=1:3
    %figure
    nurbs_plot(nurbsSurf4, subd);
    %hold on
    nurbs_plot(nurbsShape4(i),subd)
end


figure
hold on

for i=3:4
    nurbs_plot(nurbsSurf1, subd);
    nurbs_plot(nurbsShape1(i),subd)
end


for i=3:4
    nurbs_plot(nurbsSurf2, subd);
    nurbs_plot(nurbsShape2(i),subd)
end


for i=4:5
    nurbs_plot(nurbsSurf3, subd);
    nurbs_plot(nurbsShape3(i),subd)
end


for i=4:5
    nurbs_plot(nurbsSurf4, subd);
    nurbs_plot(nurbsShape4(i),subd)
end
%
elem_num = 39;
nurbsSurf5 = mesh.elem.parent(elem_num).child(1).nurbsSurf;
nurbsSurf6 = mesh.elem.parent(elem_num).child(2).nurbsSurf;
iPEN = 4:mesh.elem.parent(elem_num).nNode;

[~,iCENP,iCEN]=intersect(mesh.elem.parent(elem_num).node(iPEN),...
                        mesh.elem.parent(elem_num).child(1).node,'stable');
nurbsInd=mesh.elem.parent(elem_num).child(1).nurbsInd(iCEN);
[nurbsShape5,~]=create_nurbs_enrichment_function(nurbsSurf5,nurbsInd);

[~,iCENP,iCEN]=intersect(mesh.elem.parent(elem_num).node(iPEN),...
                        mesh.elem.parent(elem_num).child(1).node,'stable');
nurbsInd=mesh.elem.parent(elem_num).child(1).nurbsInd(iCEN);
[nurbsShape6,~]=create_nurbs_enrichment_function(nurbsSurf6,nurbsInd);
figure

hold on
nurbs_plot(nurbsSurf5,subd);
for i = 1:length(nurbsShape5)
    nurbs_plot(nurbsShape5(i),subd);
end

nurbs_plot(nurbsSurf6,subd);
for i = 1:length(nurbsShape6)
    nurbs_plot(nurbsShape6(i),subd);
end

%
nurbs_plot(nurbsSurf2,subd);
for i = 1:length(nurbsShape2)
    nurbs_plot(nurbsShape2(i),subd);
end

nurbs_plot(nurbsSurf3,subd);
for i = 1:length(nurbsShape3)
    nurbs_plot(nurbsShape3(i),subd);
end
%
%{
elem_num = 37;
nurbsSurf7 = mesh.elem.parent(elem_num).child(1).nurbsSurf;
nurbsSurf8 = mesh.elem.parent(elem_num).child(2).nurbsSurf;
iPEN = 4:mesh.elem.parent(elem_num).nNode;

elem_num=4;
nurbsSurf9=mesh.elem.parent(elem_num).child(1).nurbsSurf;
nurbsSurf10=mesh.elem.parent(elem_num).child(2).nurbsSurf;
iPEN=4:mesh.elem.parent(elem_num).nNode;
%}