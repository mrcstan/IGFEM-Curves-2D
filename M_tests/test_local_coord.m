close all
%
Xel = [0,0;
       1,0.2;
       1.2,1.2;
       0.1,0.9]';
%
%{
Xel = [-1,-1;
        1,-1;
        1, 1;
       -1, 1]';
%}
Xglo = [0.5,0.1]';
%{
xb = [-5,5];
yb = [-5,5];
Xel = rand(2,4);
Xel(1,:) = (xb(2)-xb(1))*Xel(1,:)+xb(1);
Xel(2,:) = (yb(2)-yb(1))*Xel(2,:)+yb(1);

[Xel(1,:),Xel(2,:)] = poly2ccw(Xel(1,:),Xel(2,:));

Xglo = rand(2,1);
Xglo(1) = (xb(2)-xb(1))*Xglo(1)+xb(1);
Xglo(2) = (yb(2)-yb(1))*Xglo(2)+yb(1);
while (~inpolygon(Xglo(1),Xglo(2),Xel(1,:),Xel(2,:)))
    Xglo = rand(2,1);
    Xglo(1) = (xb(2)-xb(1))*Xglo(1)+xb(1);
    Xglo(2) = (yb(2)-yb(1))*Xglo(2)+yb(1);
end
%}
patch(Xel(1,:),Xel(2,:),'r')
hold on
plot(Xglo(1),Xglo(2),'ko','markersize',10,'markerfacecolor','k')
disp('Xglo')
disp(Xglo)
Xloc = local_coord(Xglo,Xel,2);
disp('Xloc')
disp(Xloc)
N = shape_funct(Xloc,2);
Xglo2 = Xel*N;
disp('Xglo2')
disp(Xglo2)