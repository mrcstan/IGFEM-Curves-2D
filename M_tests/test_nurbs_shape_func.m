%% Using matlab quad2d function for integration of nurbs regions
% Gauss - Quadrature for NURBS element

% Preliminary settings
close all
% clear

format long e;

% Adding required toolbox pathes and reading input file
path(path, '../NURBS/nurbs_toolbox') % to add nurbs toolbox path

ftsz = 24;
subd=[30,30];
specs.colormap='default';

%
elem_num = 122;
nChilds = 4;
child_num = 1;

for i = 1:nChilds
    nurbsSurf(i)=mesh.elem.parent(elem_num).child(i).nurbsSurf;
end
iPEN=4:mesh.elem.parent(elem_num).nNode;

%
figure
hold on
for i = 1:nChilds
    specs.C = i;
    nurbs_plot(nurbsSurf(i), subd, specs);
    XX = mesh.node.coords(mesh.elem.parent(elem_num).child(i).node,:);
    for j = 1:mesh.elem.parent(elem_num).child(i).nNode
        text(XX(j,1),XX(j,2),num2str(mesh.elem.parent(elem_num).child(i).node(j)),'color','r','fontsize',16)
    end
    XX = mean(XX(1:3,:),1);
    text(XX(1),XX(2),num2str(i),'color','k','fontsize',ftsz)
end
%
xlabel('x','fontsize',ftsz)
ylabel('y','fontsize',ftsz)
zlabel('z','fontsize',ftsz)

for child_num = 1:nChilds
    [enNodes,iCENP,iCEN]=intersect(mesh.elem.parent(elem_num).node(iPEN),...
                            mesh.elem.parent(elem_num).child(child_num).node,'stable');
    nurbsInd=mesh.elem.parent(elem_num).child(child_num).nurbsInd(iCEN);
    nShape=length(nurbsInd);
    [nurbsShape,~]=create_nurbs_enrichment_function(nurbsSurf(child_num),nurbsInd);

    for i=1:nShape
        %figure
        %nurbs_plot(nurbsSurf(child_num), subd);   
        %hold on
        nurbs_plot(nurbsShape(i),subd)
        xlabel('x','fontsize',ftsz)
        ylabel('y','fontsize',ftsz)
        zlabel('z','fontsize',ftsz)
        text(mesh.node.coords(enNodes(i),1),mesh.node.coords(enNodes(i),2),num2str(enNodes(i)),'color','r','fontsize',ftsz)

    end

end


elem_num = 121;
nChilds = 4;
child_num = 2;

for i = 1:nChilds
    nurbsSurf(i)=mesh.elem.parent(elem_num).child(i).nurbsSurf;
end
iPEN=4:mesh.elem.parent(elem_num).nNode;

for i = 1:nChilds
    specs.C = i;
    nurbs_plot(nurbsSurf(i), subd, specs);
    XX = mesh.node.coords(mesh.elem.parent(elem_num).child(i).node,:);
    for j = 1:mesh.elem.parent(elem_num).child(i).nNode
        text(XX(j,1),XX(j,2),num2str(mesh.elem.parent(elem_num).child(i).node(j)),'color','r','fontsize',16)
    end
    XX = mean(XX(1:3,:),1);
    text(XX(1),XX(2),num2str(i),'color','k','fontsize',ftsz)
end

for child_num = 1:nChilds
    [enNodes,iCENP,iCEN]=intersect(mesh.elem.parent(elem_num).node(iPEN),...
                            mesh.elem.parent(elem_num).child(child_num).node,'stable');
    nurbsInd=mesh.elem.parent(elem_num).child(child_num).nurbsInd(iCEN);
    nShape=length(nurbsInd);
    [nurbsShape,~]=create_nurbs_enrichment_function(nurbsSurf(child_num),nurbsInd);

    for i=1:nShape
        %figure
        %nurbs_plot(nurbsSurf(child_num), subd);   
        %hold on
        nurbs_plot(nurbsShape(i),subd)
        xlabel('x','fontsize',ftsz)
        ylabel('y','fontsize',ftsz)
        zlabel('z','fontsize',ftsz)
        text(mesh.node.coords(enNodes(i),1),mesh.node.coords(enNodes(i),2),num2str(enNodes(i)),'color','r','fontsize',ftsz)

    end

end