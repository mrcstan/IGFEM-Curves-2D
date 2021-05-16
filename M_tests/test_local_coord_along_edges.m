close all

Xel = [6.177124344467705e-02     7.000000000000001e-02     7.000000000000001e-02;
       3.822875655532296e-02     3.000000000000000e-02     3.464101615137755e-02];
Xglo = [6.351018431618273e-02;
       3.747057779735044e-02];

edge = 3;

if (size(Xel,2) == 3)
    shape = 1;
else
    shape = 2;
end
patch(Xel(1,:),Xel(2,:),'r')
hold on
for i = 1:size(Xel,2)
    text(Xel(1,i),Xel(2,i),num2str(i),'fontsize',20)
end
plot(Xglo(1),Xglo(2),'ko','markersize',10,'markerfacecolor','k')
disp('Xglo')
disp(Xglo)
Xloc = local_coord_along_edge(Xglo,Xel,shape,edge);
disp('Xloc')
disp(Xloc)
N = shape_funct(Xloc,shape);
Xglo2 = Xel*N;
disp('Xglo2')
disp(Xglo2)