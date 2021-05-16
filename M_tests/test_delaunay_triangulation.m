close all
clear all
path(path,'/media/sf_Google_Drive/MatlabUsefulFunctions/geom_toolbox')
path(path,'/media/sf_Google_Drive/vflib/mex')
%{
XX = [4.750000000000000e-02     1.300000000000000e-01;
     4.750000000000000e-02     1.325000000000000e-01;
     4.500000000000000e-02     1.300000000000000e-01;
     4.548655596710403e-02     1.304865559671040e-01;
     4.750000000000001e-02     1.310908383759717e-01;
     4.750000000000000e-02     1.313690430684176e-01;
     4.611639864340841e-02     1.311163986434084e-01;
     4.623290000000000e-02     1.310270000000000e-01;
     4.682100000000000e-02     1.311700000000000e-01];
%}
ftsz = 30;
desiredNedges = 10;
rng('shuffle')
nNodesl2 = 6;
XX = [0,0.19;[0.15*rand(nNodesl2,1),0.2*rand(nNodesl2,1)];0.15,0.01];
constraint = [];

DT = delaunayTriangulation(XX);
E = edges(DT);
FBtri = freeBoundary(DT);
trimesh(DT.ConnectivityList,XX(:,1),XX(:,2),'color','k');
hold on
plot(XX(:,1),XX(:,2),'ro','markerfacecolor','r')
for i = 1:size(XX,1)
    text(XX(i,1),XX(i,2),num2str(i),'fontsize',ftsz);
end
axis image
set(gca,'fontsize',30)

[~,G] = adjacency_matrix(DT.ConnectivityList,nNodesl2+2);
G(logical(G)) = 1;

[Tree, pred] = graphminspantree(G,1,'Method','Kruskal');
[S, C] = graphconncomp(G,'Directed',false,'Weak',true);