%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Marcus Tan on 11/11/2013
%%% Last modified date: 11/11/2013
%%% Copyright 2013 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this function converts a NURBS curve to IGES file 
clear
path(path, '../../NURBS/nurbs_toolbox')
xi = 0.0;
xf = 0.1;
yi = 0.0;
yf = 0.05;
sumx = xi+xf;
sumy = yi+yf;
%{
yo = 0.6*sumy;
x1 = 0.2*sumx;
y1 = -0.41*sumy;
x2 = 0.4*sumx;
y2 = 1.41*sumy;
x3 = 0.6*sumx;
y3 = -0.41*sumy;
x4 = 0.8*sumx;
y4 = 1.41*sumy;
y5 = 0.4*sumy;
itrface.pts=[xi,yo; xf,y5];
itrface.contvty=[1,2];
itrface.mcf=6.0;
itrface.pt_temp = [1,20];


p=[xi,yo,0,1;
   x1,y1,0,1;
   x2,y2,0,1;
   x3,y3,0,1;
   x4,y4,0,1;
   xf,y5,0,1]';
knot=[0,0,0,1,2,3,4,4,4];
%}
ro = 0.02;
xc = 0.5*sumx;
yc = 0.0*sumy;
y1 = yc+ro;
x1 = xc-ro;
x2 = xc+ro;

w2=1/sqrt(2);

p=[x2,yc,0,1;
   x2,y1,0,w2;
   xc,y1,0,1;
   x1,y1,0,w2;
   x1,yc,0,1]';
 %p(1,:)=p(1,:).*p(4,:);
 %p(2,:)=p(2,:).*p(4,:);
 %p(3,:)=p(3,:).*p(4,:);
knot = [0 0 0 1 1 2 2 2];
itrface.nurbs=nrbmak(p,knot);

igesout({itrface.nurbs},'semicircle_ro_p02')