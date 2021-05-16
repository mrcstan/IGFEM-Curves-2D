%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Marcus Tan on 9/11/2013
%%% Last modified date: 7/2/2014
%%% Copyright 2013 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT: xi,xf,yi,yf: rectangular domain xi <= x <= xf, yi <= y <= yf
function varargout = create_channels(xi,xf,yi,yf)
% boundaries of domain
bodySource = [];
analytical = [];
u = [];
uL2 = []; % L2 norm of exact solution
% analytical body force
%{
bs=500.0;
bodySource=@(x,y) bs;
L=xf-xi;
% analytical solution
u=@(x,y) 0.5*bs*x.*(L-x);
analytical =@(x,y) [u(x,y),0.5*bs*(L-2*x),zeros(length(x),1)];
yo=0.55*(yi+yf);
channels.pts=[xi,yo; xf,yo];
channels.contvty=[1,2];
channels.mcf=0.0;
p=[xi,yo,0,1;
   xf,yo,0,1]';
knot=[0,0,1,1];
channels.nurbs=nrbmak(p,knot);
%}
% linear soln
%{
q1 = 2000;
kap = 1.0;
T2 = 20.0;
%bodySource=@(x,y) 0.0;
% analytical solution
factor = q1/kap;
const = factor*yf+T2;
u=@(x,y) -factor*y+const;
analytical =@(x,y) [u(x,y),zeros(length(x),1),-factor*ones(length(x),1)];
%{
yo=0.55*(yi+yf);
channels.pts=[xi,yo; xf,yo];
channels.contvty=[1,2];
%channels.pt_temp = [1,T2];
channels.mcf=0.0;
p=[xi,yo,0,1;
   xf,yo,0,1]';
%}
yo=0.5*(xi+xf);
channels.pts=[xi,yo; xf,yo];
channels.contvty=[1,2];
channels.mcf=0.0;
%channels.pt_temp = [1,T2];
p=[xi,yo,0,1;
   xf,yo,0,1]';

knot=[0,0,1,1];
channels.nurbs=nrbmak(p,knot);
%}
% Initialize interface.
% CASE 1: Vertical straight channel
%{
xo=0.5*(xi+xf);
channels.pts=[xo,yi; xo,yf];
channels.contvty=[1,2];
channels.mcf=1.0;
%channels.pt_temp = [1,1];
p=[xo,yi,0,1;
   xo,yf,0,1]';
knot=[0,0,1,1];
channels.nurbs=nrbmak(p,knot);
L=xf-xi;
lam=L/(channels.mcf*xo*(L-xo));
A2 = 1./(xo-L);
A1 = 1./xo;
u=@(x,y) (A1*x.*(x<=xo)+A2*(x-L).*(x>xo)).*exp(-lam*y);
ux=@(x,y) (A1*(x<=xo)+A2*(x>xo)).*exp(-lam*y);
uy=@(x,y) -lam*(A1*x.*(x<=xo)+A2*(x-L).*(x>xo)).*exp(-lam*y);
analytical=@(x,y) [u(x,y),ux(x,y),uy(x,y)];
bodySource=@(x,y) -lam^2*u(x,y);
%}
% CASE 2 Cross Channel
%{
xo=0.05;
yo=0.05;
channels.pts =[xo,yi;xo,yo;xi,yo;xf,yo;xo,yf];
channels.contvty =[1,2;2,3;2,4;2,5];
channels.pt_temp = [1,1]; % Note need to specify inlet temperature when it's not at an original node
                         % When at original node, do not need to specify or
                         % must specify correct value
% a2=a3 and mass conservation
%{
a1 = 1.0;
a2 = 0.5*(a1-4./a1);
a3 = a2;
a4 = 4./a1;
lam2 = 10;
lam1 = lam2/a1*(1+a2/a3);
lam4 = lam2/a4*(1+a2/a3);
lam3 = lam2*a2/a3;
C = exp(lam4*(yf-yo)+lam1*yo);
%}
% only mass conservation
%
a1 = 10.0;
a4 = 0.4;
lam2 = 20;

rt = roots([1,a4-a1,(a4-a1)^2/(a1*a4)]);
a2 = rt(1);
a3 = (a4-a1)^2/(a1*a4*a2);


lam1 = lam2/a1*(1+a2/a3);
lam4 = lam2/a4*(1+a2/a3);
lam3 = lam2*a2/a3;

C = exp(-lam2*xo); 

channels.mcf(1) = a1;
channels.mcf(2) = a2;
channels.mcf(3) = a3;
channels.mcf(4) = a4;
channels.lam1 = lam1;
channels.lam2 = lam2;
channels.lam3 = lam3;
channels.lam4 = lam4;
%channels.beta = beta;

%{
u=@(x,y) C*((exp(lam2*x(:)).*(x(:)<=xo)+exp(lam2*(2*xo-x(:))).*(x(:)>xo))...
        .*(exp(-lam1*y(:)).*(y(:)<=yo)+exp(lam4*(yo-y(:))-lam1*yo).*(y(:)>yo)));
ux=@(x,y) C*((lam2*exp(lam2*x(:)).*(x(:)<=xo)-lam2*exp(lam2*(2*xo-x(:))).*(x(:)>xo))...
        .*(exp(-lam1*y(:)).*(y(:)<=yo)+exp(lam4*(yo-y(:))-lam1*yo).*(y(:)>yo)));    
uy=@(x,y) C*((exp(lam2*x(:)).*(x(:)<=xo)+exp(lam2*(2*xo-x(:))).*(x(:)>xo))...
        .*(-lam1*exp(-lam1*y(:)).*(y(:)<=yo)-lam4*exp(lam4*(yo-y(:))-lam1*yo).*(y(:)>yo)));
bodySource...
  =@(x,y) C*(-lam2^2*(exp(lam2*x).*(x<=xo)+exp(lam2*(2*xo-x)).*(x>xo))...
        .*(exp(-lam1*y).*(y<=yo)+exp(lam4*(yo-y)-lam1*yo).*(y>yo))...
          -(exp(lam2*x).*(x<=xo)+exp(lam2*(2*xo-x)).*(x>xo))...
        .*(lam1^2*exp(-lam1*y).*(y<=yo)+lam4^2*exp(lam4*(yo-y)-lam1*yo).*(y>yo)));:w
%}
u=@(x,y) C*((exp(lam2*x(:)).*(x(:)<=xo)+exp(lam2*xo+lam3*(xo-x(:))).*(x(:)>xo))...
        .*(exp(-lam1*y(:)).*(y(:)<=yo)+exp(lam4*(yo-y(:))-lam1*yo).*(y(:)>yo)));
ux=@(x,y) C*((lam2*exp(lam2*x(:)).*(x(:)<=xo)-lam3*exp(lam2*xo+lam3*(xo-x(:))).*(x(:)>xo))...
        .*(exp(-lam1*y(:)).*(y(:)<=yo)+exp(lam4*(yo-y(:))-lam1*yo).*(y(:)>yo)));    
uy=@(x,y) C*((exp(lam2*x(:)).*(x(:)<=xo)+exp(lam2*xo+lam3*(xo-x(:))).*(x(:)>xo))...
        .*(-lam1*exp(-lam1*y(:)).*(y(:)<=yo)-lam4*exp(lam4*(yo-y(:))-lam1*yo).*(y(:)>yo)));
analytical=@(x,y) [u(x,y),ux(x,y),uy(x,y)];
bodySource...
  =@(x,y) C*((-lam2^2*exp(lam2*x).*(x<=xo)-lam3^2*exp(lam2*xo+lam3*(xo-x(:))).*(x>xo))...
        .*(exp(-lam1*y).*(y<=yo)+exp(lam4*(yo-y)-lam1*yo).*(y>yo))...
          -(exp(lam2*x).*(x<=xo)+exp(lam2*xo+lam3*(xo-x(:))).*(x>xo))...
        .*(lam1^2*exp(-lam1*y).*(y<=yo)+lam4^2*exp(lam4*(yo-y)-lam1*yo).*(y>yo)));
%
knot=[0,0,1,1];
p=[xo,yi;
   %xo,0.7*yi+0.3*yo;   
   %xo,0.3*yi+0.7*yo; 
   %xo,0.5*(yi+yo);
   xo,yo]';
channels.nurbs(1)=nrbmak(p,knot);
p=[xo,yo;
   %0.7*xo+0.3*xi,yo;
   %0.3*xo+0.7*xi,yo;
   %0.5*(xo+xi),yo;
   xi,yo]';
channels.nurbs(2)=nrbmak(p,knot);
p=[xo,yo;
   %0.3*xf+0.7*xo,yo;
   %0.7*xf+0.3*xo,yo;
   %0.5*(xf+xo),yo;
   xf,yo]';
channels.nurbs(3)=nrbmak(p,knot);
p=[xo,yo;
   %xo,0.3*yf+0.7*yo;
   %xo,0.7*yf+0.3*yo;
   %xo,0.5*(yf+yo);
   xo,yf]';
channels.nurbs(4)=nrbmak(p,knot);
%}
% CASE 3: 3+3+3 network
%{
sumx = xi+xf;
sumy = yi+yf;
yb1 = 0.25*sumy;
yb2 = 0.5*sumy;
yb3 = 0.75*sumy;
xj1 = 0.3*sumx;
yj1 = 0.5*sumy;
xj2 = 0.7*sumx;
yj2 = 0.5*sumy;
channels.pts = [xi,yb1;xi,yb2;xi,yb3;xj1,yj1;xj2,yj2;xf,yb1;xf,yb2;xf,yb3];
channels.contvty = [1,4;2,4;3,4;4,5;4,5;4,5;5,6;5,7;5,8];
channels.mcf = [1.0,1.0,1.0,  0.7,1.4,0.7, 0.5,1.5,0.5];

knot=[0,0,0,1,1,1];
% left sections
x1 = 0.2*sumx;
y1 = 0.35*sumy;
p=[xi,yb1;
   x1,y1;
   xj1,yj1]';
channels.nurbs(1)=nrbmak(p,knot);
y1 = 0.6*sumy;
p=[xi,yb2;
   x1,y1;
   xj1,yj1]';
channels.nurbs(2)=nrbmak(p,knot);
y1 = 0.65*sumy;
p=[xi,yb3;
   x1,y1;
   xj1,yj1]';
channels.nurbs(3)=nrbmak(p,knot);
% mid sections
x1 = 0.5*sumx;
y1 = 0.1*sumy;
p=[xj1,yj1;
   x1,y1;
   xj2,yj2]';
channels.nurbs(4)=nrbmak(p,knot);
x1 = 0.5*sumx;
y1 = 0.45*sumy;
p=[xj1,yj1;
   x1,y1;
   xj2,yj2]';
channels.nurbs(5)=nrbmak(p,knot);
x1 = 0.5*sumx;
y1 = 0.95*sumy;
p=[xj1,yj1;
   x1,y1;
   xj2,yj2]';
channels.nurbs(6)=nrbmak(p,knot);
% right sections
x1 = 0.85*sumx;
y1 = 0.3*sumy;
p=[xj2,yj2;
   x1,y1;
   xf,yb1]';
channels.nurbs(7)=nrbmak(p,knot);
channels.nNurbs=length(channels.nurbs);
x1 = 0.85*sumx;
y1 = 0.4*sumy;
p=[xj2,yj2;
   x1,y1;
   xf,yb2]';
channels.nurbs(8)=nrbmak(p,knot);
x1 = 0.85*sumx;
y1 = 0.8*sumy;
p=[xj2,yj2;
   x1,y1;
   xf,yb3]';
channels.nurbs(9)=nrbmak(p,knot);
%}
% CASE 4: 2+3+1 network
%{
sumx = xi+xf;
sumy = yi+yf;
yb1 = 0.2*sumy;
yb2 = 0.5*sumy;
yb3 = 0.8*sumy;
xj1 = 0.29*sumx;
yj1 = 0.5*sumy;
xj2 = 0.69*sumx;
yj2 = 0.5*sumy;
channels.pts = [xi,yb1;xi,yb3;xj1,yj1;xj2,yj2;xf,yb2];
channels.contvty = [1,3;2,3;3,4;3,4;3,4;4,5];
channels.mcf = [1.5,1.5, 1,0,1.0,1.0, 3.0]*2.0;
channels.pt_temp = [1,20;2,20];

knot=[0,0,0,1,1,1];
% left sections
x1 = 0.2*sumx;
y1 = 0.35*sumy;
p=[xi,yb1;
   x1,y1;
   xj1,yj1]';
channels.nurbs(1)=nrbmak(p,knot);
y1 = 0.65*sumy;
p=[xi,yb3;
   x1,y1;
   xj1,yj1]';
channels.nurbs(2)=nrbmak(p,knot);
% mid sections
x1 = 0.5*sumx;
y1 = 0.1*sumy;
p=[xj1,yj1;
   x1,y1;
   xj2,yj2]';
channels.nurbs(3)=nrbmak(p,knot);
x1 = 0.5*sumx;
y1 = 0.45*sumy;
p=[xj1,yj1;
   x1,y1;
   xj2,yj2]';
channels.nurbs(4)=nrbmak(p,knot);
x1 = 0.5*sumx;
y1 = 0.95*sumy;
p=[xj1,yj1;
   x1,y1;
   xj2,yj2]';
channels.nurbs(5)=nrbmak(p,knot);
% right section
x1 = 0.85*sumx;
y1 = 0.4*sumy;
p=[xj2,yj2;
   x1,y1;
   xf,yb2]';
channels.nurbs(6)=nrbmak(p,knot);
%}
% CASE 5: wavy channel
%{
sumx = xi+xf;
sumy = yi+yf;
yo = 0.6*sumy;
x1 = 0.2*sumx;
y1 = 0.09*sumy;
x2 = 0.4*sumx;
y2 = 0.91*sumy;
x3 = 0.6*sumx;
y3 = 0.09*sumy;
x4 = 0.8*sumx;
y4 = 0.91*sumy;
y5 = 0.4*sumy;
channels.pts=[xi,yo; xf,y5];
channels.contvty=[1,2];
channels.mcf = 6.0;
channels.pt_temp = [1,20.0];
channels.kapf = 1.0;
channels.Tin = 20.0;

p=[xi,yo,0,1;
   x1,y1,0,1;
   x2,y2,0,1;
   x3,y3,0,1;
   x4,y4,0,1;
   xf,y5,0,1]';
knot=[0,0,0,1,2,3,4,4,4];
channels.nurbs=nrbmak(p,knot);

% for conforming mesh
% channels.region_region_seg = sparse(1,2,1); 
channels.seg_region_region = [1,2]; 
%}
% CASE 6: wavy channel with high_curvature
%{
sumx = xi+xf;
sumy = yi+yf;
yo = 0.6*sumy;
x1 = 0.1*sumx;
y1 = -0.45*sumy;
x2 = 0.25*sumx;
y2 = 1.45*sumy;
x3 = 0.35*sumx;
y3 = -0.45*sumy;
x4 = 0.5*sumx;
y4 = 1.45*sumy;
x5 = 0.65*sumx;
y5 = -0.45*sumy;
x6 = 0.75*sumx;
y6 = 1.45*sumy;
x7 = 0.9*sumx;
y7 = -0.4*sumy;
y8 = 0.4*sumy;
channels.pts=[xi,yo; xf,y5];
channels.contvty=[1,2];
channels.mcf = 30.0;
channels.pt_temp = [1,20.0];
channels.kapf = 1.0;
channels.Tin = 20.0;

p=[xi,yo;
   x1,y1;
   x2,y2;
   x3,y3;
   x4,y4;
   x5,y5;
   x6,y6;
   x7,y7;
   xf,y8]';
knot=[0,0,0,1,2,3,4,5,6,7,7,7];
channels.nurbs = nrbmak(p,knot);

% for conforming mesh
channels.region_region_seg = sparse(1,2,1); 
%}
% CASE 7: wavy channel with greater amplitude
%{
sumx = xi+xf;
sumy = yi+yf;
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
channels.pts=[xi,yo; xf,y5];
channels.contvty=[1,2];
channels.mcf = 18.0;
channels.pt_temp = [1,20];
channels.region_region_seg = sparse(1,2,1);

p=[xi,yo,0,1;
   x1,y1,0,1;
   x2,y2,0,1;
   x3,y3,0,1;
   x4,y4,0,1;
   xf,y5,0,1]';
knot=[0,0,0,1,2,3,4,4,4];
channels.nurbs=nrbmak(p,knot);
%}
% CASE 8: wavy channel with more oscillations
%{
sumx = xi+xf;
sumy = yi+yf;
yo = 0.6*sumy;
x1 = 0.1*sumx;
y1 = -0.4*sumy;
x2 = 0.2*sumx;
y2 = 1.4*sumy;
x3 = 0.3*sumx;
y3 = -0.4*sumy;
x4 = 0.4*sumx;
y4 = 1.4*sumy;
x5 = 0.5*sumx;
y5 = -0.4*sumy;
x6 = 0.6*sumx;
y6 = 1.4*sumy;
x7 = 0.7*sumx;
y7 = -0.4*sumy;
x8 = 0.8*sumx;
y8 = 1.4*sumy;
x9 = 0.9*sumx;
y9 = -0.4*sumy;
y10 = 0.6*sumy;
channels.pts=[xi,yo; xf,y10];
channels.contvty=[1,2];
channels.mcf=10.0;
channels.pt_temp = [1,20];


p=[xi,yo;
   x1,y1;
   x2,y2;
   x3,y3;
   x4,y4;
   x5,y5;
   x6,y6;
   x7,y7;
   x8,y8;
   x9,y9;
   xf,y10]';
knot=[0,0,0,1,2,3,4,5,6,7,8,9,9,9];
channels.nurbs = nrbmak(p,knot);

%}
%{
sumx = xi+xf;
sumy = yi+yf;
yo = 0.6*sumy;
x1 = 0.2*sumx;
y1 = 0.09*sumy;
x2 = 0.4*sumx;
y2 = 0.91*sumy;
x3 = 0.6*sumx;
y3 = 0.09*sumy;
x4 = 0.8*sumx;
y4 = 0.91*sumy;
y5 = 0.4*sumy;
channels.pts=[xi,yo; xf,y5];
channels.contvty=[1,2];
channels.mcf=500.0;
channels.pt_temp = [1,20];

p=[xi,yo,0,1;
   x1,y1,0,1;
   x2,y2,0,1;
   x3,y3,0,1;
   x4,y4,0,1;
   xf,y5,0,1]';
knot=[0,0,0,1,2,3,4,4,4];
bodySource = @(x,y) 5e5;
%}
% CASE 9: wavy network
%{
sumx = xi+xf;
sumy = yi+yf;
yb1 = 0.6*sumy;
yb2 = 0.4*sumy;
xj1 = 0.11*sumx;
yj1 = 0.51*sumy;
xj2 = 0.89*sumx;
yj2 = 0.51*sumy;
xm1 = 0.3*sumx;
xm2 = 0.5*sumx;
xm3 = 0.7*sumx;
ym11 = 1.0*sumy;
ym12 = 0.65*sumy;
ym13 = 1.0*sumy;
ym21 = 0.75*sumy;
ym22 = 0.5*sumy;
ym23 = 0.75*sumy;
ym31 = 0.25*sumy;
ym32 = 0.5*sumy;
ym33 = 0.25*sumy;
ym41 = 0.0*sumy;
ym42 = 0.35*sumy;
ym43 = 0.0*sumy;
channels.pts = [xi,yb1;xj1,yj1;xj2,yj2;xf,yb2];
channels.contvty = [1,2;2,3;2,3;2,3;2,3;3,4];
                                       
channels.mcf = [6.0 ,1.5, 1.5, 1.5, 1.5, 6.0];
channels.pt_temp = [1,20];
knot=[0,0,1,1];
% left sections
p=[xi,yb1;
   xj1,yj1]';
channels.nurbs(1)=nrbmak(p,knot);
% mid sections
knot = [0,0,0,1,2,3,3,3];
p=[xj1,yj1;
   xm1,ym11;
   xm2,ym12;
   xm3,ym13;
   xj2,yj2]';
channels.nurbs(2)=nrbmak(p,knot);
p=[xj1,yj1;
   xm1,ym21;
   xm2,ym22;
   xm3,ym23;
   xj2,yj2]';
channels.nurbs(3)=nrbmak(p,knot);
p=[xj1,yj1;
   xm1,ym31;
   xm2,ym32;
   xm3,ym33;
   xj2,yj2]';
channels.nurbs(4)=nrbmak(p,knot);
p=[xj1,yj1;
   xm1,ym41;
   xm2,ym42;
   xm3,ym43;
   xj2,yj2]';
channels.nurbs(5)=nrbmak(p,knot);
% right section
knot = [0,0,1,1];
p=[xj2,yj2;
   xf,yb2]';
channels.nurbs(6)=nrbmak(p,knot);
% for conforming mesh
channels.seg_region_region = [1,5;1,2;2,3;3,4;4,5;1,5];
%}
%{
diam = 0.0005;
nu = 1e-6;
            
heatCapacity = 3494; % in SI units
[pressure,mass] = network_pressure_mass_flow_rate(channels.contvty,...
                                                           channels.nurbs,...
                                                           diam,...
                                                           [],...
                                                           nu,...
                                                           1,...
                                                           6/heatCapacity,...
                                                           [],...
                                                           4,...
                                                           0,...
                                                           'circular');
plot_channel_network(channels.contvty,channels.nurbs,pressure,mass)

channels.mcf = mass*heatCapacity;
%}
% CASE 10: serpentine channel
%{
sumx = xi+xf;
sumy = yi+yf;
ys = 0.1*sumy;
x1 = 0.89*sumx;
y1 = 0.039*sumy;
x2 = 0.89*sumx;
y2 = 0.69*sumy;
x3 = 0.3*sumx;
%x3 = 0.45*sumx;
y3 = 0.69*sumy;
x4 = 0.3*sumx;
y4 = 0.5*sumy;
x5 = 0.7*sumx;
y5 = 0.6*sumy;
x6 = 0.6*sumx;
y6 = 0.25*sumy;
x7 = 0.19*sumx;
y7 = 0.25*sumy;
x8 = 0.19*sumx;
y8 = 0.98*sumx;
ye = 0.9*sumy;
channels.pts=[xi,ys; xf,ye];
channels.contvty=[1,2];
%channels.mcf=1000.0;
%channels.pt_temp = [1,20];
channels.mcf = 2.5;
channels.pt_temp = [1,0];
p=[xi,ys;
   x1,y1;
   x2,y2;
   x3,y3;
   x4,y4;
   x5,y5;
   x6,y6;
   x7,y7;
   x8,y8;
   xf,ye]';
knot=[0,0,0,1,2,3,4,5,6,7,8,8,8];
%bodySource = @(x,y) 5e5;

channels.nurbs=nrbmak(p,knot);
%}
% CASE 11: curvature effect on solution (NOT WORKING !!!!)
%{
sumx = xi+xf;
sumy = yi+yf;
yo = 0.5*sumy;
xm = 0.5*sumx;
ym = 0.8*sumy;
channels.pts=[xi,yo; xf,yo];
channels.contvty=[1,2];
channels.mcf=10.0;
channels.pt_temp = [1,20];


p=[xi,yo,0,1;
   xm,ym,0,1;
   xf,yo,0,1]';
knot=[0,0,0,1,1,1];
channels.nurbs=nrbmak(p,knot);
%}
% CASE 12: Semicircular channel CCW flow
%
n=3;
ro = 0.4*(xf-xi);
sumx = xi+xf;
sumy = yi+yf;
xc = 0.5*sumx;
yc = 0.0*sumy;
y1 = yc+ro;
x1 = xc-ro;
x2 = xc+ro;
channels.pts=[xf,yc; xi,yc];
channels.contvty=[1,2];
channels.mcf = 90.0;
channels.xc = xc;
channels.yc = yc;
channels.ro = ro;
channels.pt_temp = [1,1];
channels.region_region_seg = sparse(1,2,1);

w2=1/sqrt(2);

p=[x2,yc,0,1;
   x2,y1,0,w2;
   xc,y1,0,1;
   x1,y1,0,w2;
   x1,yc,0,1]';
 p(1,:)=p(1,:).*p(4,:);
 p(2,:)=p(2,:).*p(4,:);
 p(3,:)=p(3,:).*p(4,:);
knot = [0 0 0 1 1 2 2 2];
channels.nurbs=nrbmak(p,knot);

c1 = 1/ro^n;
lam = 2*n/channels.mcf;
c2 = c1*ro^(2*n);
%
u = @(x,y) u_func(x,y,xc,yc,ro,c1,c2,n,lam);
ux = @(x,y) ux_func(x,y,xc,yc,ro,c1,c2,n,lam);
uy = @(x,y) uy_func(x,y,xc,yc,ro,c1,c2,n,lam);
%ur = @(x,y) ur_func(x,y,xc,yc,ro,c1,c2,m,lam);
analytical=@(x,y) [u(x,y),ux(x,y),uy(x,y)];
bodySource=@(x,y) source_func(x,y,xc,yc,ro,c1,c2,n,lam);
uL2 = u_exact_L2_norm(ro,c1,c2,n,lam,0.5*(xf-xi));
%
% CASE 13: Semicircular channel CW flow (NOT WORKING !!!!)
%{
n=3;
ro = 0.041;
sumx = xi+xf;
sumy = yi+yf;
xc = 0.5*sumx;
yc = 0.0*sumy;
y1 = yc+ro;
x1 = xc-ro;
x2 = xc+ro;
channels.pts=[xi,yc; xf,yc];
channels.contvty=[1,2];
channels.mcf=10;
%channels.pt_temp = [1,20];
channels.region_region_seg = sparse(1,2,1);

w2=1/sqrt(2);

p=[x1,yc,0,1;
   x1,y1,0,w2;
   xc,y1,0,1;
   x2,y1,0,w2;
   x2,yc,0,1]';
 p(1,:)=p(1,:).*p(4,:);
 p(2,:)=p(2,:).*p(4,:);
 p(3,:)=p(3,:).*p(4,:);
knot = [0 0 0 1 1 2 2 2];
channels.nurbs=nrbmak(p,knot);
c1 = 1/ro^n;
lam = 2*n/channels.mcf;
c2 = c1*ro^(2*n);
%
u = @(x,y) u_func(x,y,xc,yc,ro,c1,c2,n,lam);
ux = @(x,y) ux_func(x,y,xc,yc,ro,c1,c2,n,lam);
uy = @(x,y) uy_func(x,y,xc,yc,ro,c1,c2,n,lam);
%ur = @(x,y) ur_func(x,y,xc,yc,ro,c1,c2,m,lam);
analytical=@(x,y) [u(x,y),ux(x,y),uy(x,y)];
bodySource=@(x,y) source_func(x,y,xc,yc,ro,c1,c2,n,lam);
%}
% CASE 14: 6 channels - comparison with experiments
%{
%xshift = 0.0142725864425302;
xshift = 0.0;
x1 = 0.013 + xshift;
x2 = 0.0478 + xshift;
x3 = 0.0826 + xshift;
x4 = 0.1174 + xshift;
x5 = 0.1522 + xshift;
x6 = 0.187 + xshift;
channels.pts =[x1,yi;x1,yf;
              x2,yi;x2,yf;
              x3,yi;x3,yf;
              x4,yi;x4,yf;
              x5,yi;x5,yf;
              x6,yi;x6,yf];
channels.contvty =[1,2;3,4;5,6;7,8;9,10;11,12];
channels.pt_temp = [1,23.0;
                   3,23.0;
                   5,23.0;
                   7,23.0;
                   9,23.0;
                   11,23.0];

mdot = 1.747/6.0;
for i = 1:6                
    channels.mcf(i) = mdot;
    channels.kapf(i) = 0.419;
    channels.Tin(i) = 23.0;
    channels.length(i) = 0.15;
end               
               
knot = [0,0,1,1];
p = [x1,yi;
     x1,yf]';
channels.nurbs(1) = nrbmak(p,knot);

p = [x2,yi;
     x2,yf]';
channels.nurbs(2) = nrbmak(p,knot);

p = [x3,yi;
     x3,yf]';
channels.nurbs(3) = nrbmak(p,knot);

p = [x4,yi;
     x4,yf]';
channels.nurbs(4) = nrbmak(p,knot);

p = [x5,yi;
     x5,yf]';
channels.nurbs(5) = nrbmak(p,knot);

p=[x6,yi;
   x6,yf]';
channels.nurbs(6) = nrbmak(p,knot);
%bodySource = @(x,y) 500.0;
%}
%  CASE 14: Curved channel consisting of circular arcs for verification with FLUENT
%{
w2=1/sqrt(2);
channels.pts = [0,0.01;0.12,0.05];
channels.contvty = [1,2];
channels.pt_temp = [1,23.0];             
channels.mcf = 1.747/6.0;
channels.kapf = 0.419;
channels.Tin = 23.0;              

p = [0,0.01,0,1;
     0.02,0.01,0,w2;
     0.02,0.03,0,1;
     0.02,0.05,0,w2;
     0.04,0.05,0,1;
     0.06,0.05,0,w2;
     0.06,0.03,0,1;
     0.06,0.01,0,w2;
     0.08,0.01,0,1;
     0.1,0.01,0,w2;
     0.1,0.03,0,1;
     0.1,0.05,0,w2;
     0.12,0.05,0,1]';
p(1,:)=p(1,:).*p(4,:);
p(2,:)=p(2,:).*p(4,:);

knot = [0,0,0,1,1,2,2,3,3,4,4,5,5,6,6,6]; 
channels.nurbs=nrbmak(p,knot);
channels.region_region_seg = sparse(1,2,1);
%}
% CASE 5a: Problematic Channel 1
%{
channels.pts=[0,0.0300000004200000; 0.100000001000000,0.0200000002800000];
channels.contvty=[1,2];
channels.mcf=90.0;
channels.pt_temp = [1,20];

p = [0,0.0300000004200000;
     0.0200000002000000,0.0174518453438583;
     0.0400000004000000,0.0429559737029367;
     0.0600000006000000,0.0153936890764418;
     0.0800000008000000,0.0450424547849753;
     0.100000001000000,0.0200000002800000]';
knot=[0,0,0,0.250000000000000,0.500000000000000,0.750000000000000,1,1,1;];
channels.nurbs=nrbmak(p,knot);

%}
% CASE 15: Debug 3 intersections on one regular element
%{
channels.pts = [0,0.10;0.10,0];
channels.contvty = [1,2];
channels.pt_temp = [1,21.5];             
channels.mcf = 1.747;
channels.kapf = 0.419;
channels.Tin = 21.5;  
p = [0,0.1;
     0.1,0.1;
     0.1,0]';
knot = [0,0,1,2,2];
channels.nurbs = nrbmak(p,knot);
%}
% CASE 16: Spiral channel verification and validation
%{
channels.pts = [0,0.189;0.15,0.011];
channels.contvty = [1,2];
channels.pt_temp = [1,21.5];             
channels.mcf = 1.747;
channels.kapf = 0.419;
channels.Tin = 21.5;  
p = [0,0.189;
     0.141,0.189;
     0.141,0.033;
     0.0285,0.033;
     0.0285,0.144;
     0.104,0.144;
     0.104,0.078;
     0.065,0.078;
     0.065,0.1;
     0.085,0.1;
     0.085,0.122;
     0.0465,0.122;
     0.0465,0.055;
     0.1215,0.055;
     0.1215,0.166;
     0.009,0.166;
     0.01,0.011;
     0.15,0.011]';
knot = [0,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,17];
channels.nurbs = nrbmak(p,knot);
%}
%{
channels.pts =   [0,0.189;
                 0.141,0.189;
                 0.141,0.033;
                 0.0285,0.033;
                 0.0285,0.144;
                 0.104,0.144;
                 0.104,0.078;
                 0.065,0.078;
                 0.065,0.1;
                 0.085,0.1;
                 0.085,0.122;
                 0.0465,0.122;
                 0.0465,0.055;
                 0.1215,0.055;
                 0.1215,0.166;
                 0.009,0.166;
                 0.01,0.011;
                 0.15,0.011];
channels.contvty = [1,2;
                   2,3;
                   3,4;
                   4,5;
                   5,6;
                   6,7;
                   7,8;
                   8,9;
                   9,10;
                   10,11;
                   11,12;
                   12,13;
                   13,14;
                   14,15;
                   15,16;
                   16,17;
                   17,18];

channels.pt_temp = [1,21.5]; 
channels.mcf = 1.747*ones(17,1);
channels.kapf = 0.419;
channels.Tin = 21.5;   

knot = [0,0,1,1];
for i = 1:17
    channels.nurbs(i) = nrbmak([channels.pts(i,:);channels.pts(i+1,:)]',knot);
end
%}

% CASE 17: parallel channel verification and validation
%
% remember to change p whenever you change channels.pts
%
% original locations of control points
%{
channels.pts =   [0,    0.189;
                 0.01, 0.189;
                 0.01, 0.16675;
                 0.139,0.16675;
                 0.01, 0.1445;
                 0.139,0.1445;
                 0.01, 0.12225;
                 0.139,0.12225;
                 0.01, 0.1;
                 0.139,0.1;
                 0.01, 0.07775;
                 0.139,0.07775;
                 0.01, 0.0555;
                 0.139,0.0555;
                 0.01, 0.03325;
                 0.139,0.03325;
                 0.139,0.011;
                 0.15, 0.011];
%
% modified locations of control points
% for getting mesh that has the inlet coinciding with original node,
% let xi=0.0,xf=0.15,yi=0,yf=0.20, nx = 30, ny = 40
% for getting conforming mesh,
% let xi=0,xf=0.15,yi=-0.00125,yf=0.20125, nx = 15, ny = 18
%{
channels.pts =   [0,    0.19;
                 0.01, 0.19;
                 0.01, 0.1675;
                 0.14,0.1675;
                 0.01, 0.145;
                 0.14,0.145;
                 0.01, 0.1225;
                 0.14,0.1225;
                 0.01, 0.1;
                 0.14,0.1;
                 0.01, 0.0775;
                 0.14,0.0775;
                 0.01, 0.055;
                 0.14,0.055;
                 0.01, 0.0325;
                 0.14,0.0325;
                 0.14,0.01;
                 0.15, 0.01];
 %}
channels.contvty = [1,2;
                   2,4;
                   3,4;
                   5,6;
                   7,8;
                   9,10;
                   11,12;
                   13,14;
                   15,16;
                   15,17;
                   2,3;
                   3,5;
                   5,7;
                   7,9;
                   9,11;
                   11,13;
                   13,15;
                   4,6;
                   6,8;
                   8,10;
                   10,12;
                   12,14;
                   14,16;
                   16,17;
                   17,18];

channels.pt_temp = [1,21.5]; 
channels.kapf = 0.419;
channels.Tin = 21.5;   

knot = [0,0,1,1];
for i = 1:size(channels.contvty,1)
    channels.nurbs(i) = nrbmak(channels.pts(channels.contvty(i,:),:)',knot);
end
%
p = [0.01,0.189;
     0.139,0.189;
     0.139,0.16675]';
channels.nurbs(2) = nrbmak(p,[0,0,1,2,2]);
p = [0.01,0.03325;
     0.01,0.011;
     0.139,0.011]';
channels.nurbs(10) = nrbmak(p,[0,0,1,2,2]);
%
%modified locations of control points
%{
p = [0.01,0.19;
     0.14,0.19;
     0.14,0.1675]';
channels.nurbs(2) = nrbmak(p,[0,0,1,2,2]);
p = [0.01,0.0325;
     0.01,0.01;
     0.14,0.01]';
channels.nurbs(10) = nrbmak(p,[0,0,1,2,2]);
%}
%
diam = 7.5e-4; % m
nu = 3.405e-6; % viscosity
heatCapacity = 3494; % J/kg.K
massin = 5e-4; % kg/s
[pressure,mass] = network_pressure_mass_flow_rate(channels.contvty,...
                                                    channels.nurbs,...
                                                    diam,...
                                                    [],...
                                                    nu,...
                                                    1,...
                                                    massin,...
                                                    [],...
                                                    18,...
                                                    0,...
                                                    'circular');
plot_channel_network(channels.contvty,channels.nurbs,pressure,mass*1e6*60/1065)

channels.mcf = mass*heatCapacity;
%}
% CASE 18: SUPG study with cross channels
%{
xo=0.05;
yo=0.05;
channels.pts =[xo,yi;xo,yo;xi,yo;xf,yo;xo,yf];
channels.contvty =[1,2;2,3;2,4;2,5];
channels.pt_temp = [1,21.5]; 

knot=[0,0,1,1];
p=[xo,yi;
   xo,yo]';
channels.nurbs(1)=nrbmak(p,knot);
p=[xo,yo;
   xi,yo]';
channels.nurbs(2)=nrbmak(p,knot);
p=[xo,yo;
   xf,yo]';
channels.nurbs(3)=nrbmak(p,knot);
p=[xo,yo;
   xo,yf]';
channels.nurbs(4)=nrbmak(p,knot);

diam = 7.5e-4; % m
nu = 3.405e-6; % viscosity
heatCapacity = 3494; % J/kg.K
massin = 20e-4; % kg/s
[pressure,mass] = network_pressure_mass_flow_rate(channels.contvty,...
                                                    channels.nurbs,...
                                                    diam,...
                                                    [],...
                                                    nu,...
                                                    1,...
                                                    massin,...
                                                    [],...
                                                    [3,4,5],...
                                                    [0,0,0],...
                                                    'circular');
plot_channel_network(channels.contvty,channels.nurbs,pressure,mass*1e6*60/1065)

channels.mcf = mass*heatCapacity;
%}
% CASE 19: SUPG study with one channel loop
%{
x1 = 0.125 * (xi + xf);
y1 = 0.125 * (yi + yf);
x2 = 0.875 * (xi + xf);
y2 = 0.875 * (yi + yf);
channels.pts =[xi,y2;
              x1,y2;
              x2,y1;
              xf,y1];
channels.contvty =[1,2;
                  2,3;
                  2,3;
                  3,4];
channels.pt_temp = [1,21.5]; 

knot = [0,0,1,1];
p = [xi,y2;
     x1,y2]';
channels.nurbs(1)=nrbmak(p,knot);
p = [x2,y1;
     xf,y1]';
channels.nurbs(4)=nrbmak(p,knot);
knot = [0,0,1,2,2];
p = [x1,y2;
     x1,y1;
     x2,y1]';
channels.nurbs(2)=nrbmak(p,knot);
p = [x1,y2;
     x2,y2;
     x2,y1]';
channels.nurbs(3)=nrbmak(p,knot);

diam = 7.5e-4; % m
nu = 3.405e-6; % viscosity
heatCapacity = 3494; % J/kg.K
massin = 5e-4; % 5e-4 kg/s = 28 ml/min
%massin = 1.775e-5/2; % 1.775e-5 kg/s = 1 ml/m
[pressure,mass] = network_pressure_mass_flow_rate(channels.contvty,...
                                                    channels.nurbs,...
                                                    diam,...
                                                    [],...
                                                    nu,...
                                                    1,...
                                                    massin,...
                                                    [],...
                                                    4,...
                                                    0,...
                                                    'circular');
plot_channel_network(channels.contvty,channels.nurbs,pressure,mass*1e6*60/1065)

channels.mcf = mass*heatCapacity;
%}
% CASE 19: SUPG study with two channel loop
%{
x1 = 0.125 * (xi + xf);
y1 = 0.125 * (yi + yf);
y2 = 0.5 * (yi + yf);
x3 = 0.875 * (xi + xf);
y3 = 0.875 * (yi + yf);
channels.pts =[xi,y3;
              x1,y3;
              x1,y2;
              x3,y2;
              x3,y1;
              xf,y1];
channels.contvty = [1,2;
                   2,3;
                   2,4;
                   3,4;
                   3,5;
                   4,5;
                   5,6];
channels.pt_temp = [1,21.5]; 

knot = [0,0,1,1];
p = [xi,y3;
     x1,y3]';
channels.nurbs(1)=nrbmak(p,knot);
p = [x1,y3;
     x1,y2]';
channels.nurbs(2)=nrbmak(p,knot);
p = [x1,y2;
     x3,y2]';
channels.nurbs(4)=nrbmak(p,knot);
p = [x3,y2;
     x3,y1]';
channels.nurbs(6)=nrbmak(p,knot);
p = [x3,y1;
     xf,y1]';
channels.nurbs(7)=nrbmak(p,knot);

knot = [0,0,1,2,2];
p = [x1,y3;
     x3,y3;
     x3,y2]';
channels.nurbs(3)=nrbmak(p,knot);
p = [x1,y2;
     x1,y1;
     x3,y1]';
channels.nurbs(5)=nrbmak(p,knot);

diam = 7.5e-4; % m
nu = 3.405e-6; % viscosity
heatCapacity = 3494; % J/kg.K
massin = 5e-4; % 5e-4 kg/s = 28 ml/min
%massin = 1.775e-5; % 1.775e-5 kg/s = 1 ml/min
[pressure,mass] = network_pressure_mass_flow_rate(channels.contvty,...
                                                    channels.nurbs,...
                                                    diam,...
                                                    [],...
                                                    nu,...
                                                    1,...
                                                    massin,...
                                                    [],...
                                                    6,...
                                                    0,...
                                                    'circular');
plot_channel_network(channels.contvty,channels.nurbs,pressure,mass*1e6*60/1065)

channels.mcf = mass*heatCapacity;
%}
% CASE 20: SUPG study with single channel
%{
%yo = 0.5 * (yi + yf);
channels.pts =[xi,yf;
              xf,yi];
channels.contvty =[1,2];
channels.pt_temp = [1,21.5]; 
channels.Tin = 21.5;  
channels.kapf = 0.419;

knot = [0,0,1,1];

channels.nurbs=nrbmak(channels.pts',knot);


heatCapacity = 3494; % J/kg.K
massin = 50e-4; % 5e-4 kg/s = 28 ml/min
%massin = 1.775e-5; % 1.775e-5 kg/s = 1 ml/min

channels.mcf = massin*heatCapacity;
%} 
% CASE 20: Spiral channel with negative temperature
%{
p = [0,0.0999812256045221,0.0497554434077493,0.0200946908656795,0.0201606835999123,0.100000000000000;
     0.0100000000000000,0.00496289785488824,0.0488927960800881,0.0476291666915558,0.0815447844616246,0.0900000000000000];
channels.pts=[p(1,1),p(2,1); p(1,end),p(2,end)];
channels.contvty=[1,2];
heatCapacity = 3494; % J/kg.K
massin = 5e-4; % kg/s
channels.mcf = massin*heatCapacity;
channels.pt_temp = [1,21.5];
knot=[0,0,0,1,2,3,4,4,4];
%bodySource = @(x,y) 5e5;

channels.nurbs=nrbmak(p,knot);
%}
%%
channels.mcf = channels.mcf';
channels.model = 1; % model: 1 - dT/ds model
                    %        2 - h(T(s)-Tin) model
channels.nNurbs=length(channels.nurbs);

% for mex function to find intersections
channels.kind = ones(channels.nNurbs,1);
for i = 1:channels.nNurbs
    if (any(channels.nurbs(i).coefs(4,:) ~= 1))
        channels.kind(i) = 2;
        break
    end
end
channels.rows{1} = [1,2];
channels.rows{2} = [1,2,4];

% calculate nurbs arc length 
channels.length = zeros(channels.nNurbs,1);
ngpts = 10;
for i = 1:channels.nNurbs
    if (channels.nurbs(i).order == 2)
        channels.segLength{i} = sqrt(sum((channels.nurbs(i).coefs(1:2,2:end)...
                                              -channels.nurbs(i).coefs(1:2,1:end-1)).^2,1))';
        channels.cumSumSegLengths{i} = [0;cumsum(channels.segLength{i})];
        channels.length(i) = channels.cumSumSegLengths{i}(end);
        channels.lineSegs{i} = [channels.nurbs(i).coefs(1:2,1:end-1);...
                                channels.nurbs(i).coefs(1:2,2:end)]';  
    else
        channels.length(i) = nurbs_arc_length(ngpts,channels.nurbs(i));
    end
end


if (channels.model == 2)
    channels.pt_temp = [];
end
% find the junctions (terminal points where more than one line source 
% segments intersect) of the line source segments
channels.junc = branching_points(channels.contvty);                                 
channels.kinks = kinks(channels.nurbs);
varargout{1} = channels;
varargout{2} = bodySource;
varargout{3} = analytical;
varargout{4} = u;
varargout{5} = uL2;
end

function u = u_func(x,y,xc,yc,ro,c1,c2,n,lam)
    [t,r] = cart2pol(x-xc,y-yc); 
    len = length(t(:));
    u = zeros(len,1);
    for i = 1:len
        if (r(i) <= ro)
            u(i) = c1*r(i).^n.*exp(-lam*t(i));
        else
            u(i) = c2*r(i).^(-n).*exp(-lam*t(i));
        end
    end
end

function ux = ux_func(x,y,xc,yc,ro,c1,c2,n,lam)
    [t,r] = cart2pol(x-xc,y-yc);
    len = length(t(:));
    ux = zeros(len,1);
    for i = 1:len
        if (r(i) < ro)
            ux(i) = c1*r(i).^(n-1).*exp(-lam*t(i)).*(n*cos(t(i))+lam*sin(t(i)));
        else
            ux(i) = c2*r(i).^(-n-1).*exp(-lam*t(i)).*(-n*cos(t(i))+lam*sin(t(i)));
        end
    end
end

function uy = uy_func(x,y,xc,yc,ro,c1,c2,n,lam)
    [t,r] = cart2pol(x-xc,y-yc);
    len = length(t(:));
    uy = zeros(len,1);
    for i = 1:len
        if (r(i) < ro)
            uy(i) = c1*r(i).^(n-1).*exp(-lam*t(i)).*(n*sin(t(i))-lam*cos(t(i)));
        else
            uy(i) = c2*r(i).^(-n-1).*exp(-lam*t(i)).*(-n*sin(t(i))-lam*cos(t(i)));
        end
    end
end

function bodySource = source_func(x,y,xc,yc,ro,c1,c2,n,lam)
    [t,r] = cart2pol(x-xc,y-yc);
    len = length(t);
    bodySource = zeros(len,1);
    for i = 1:len
        if (r(i) < ro)
            bodySource(i) = -(n^2+lam^2)*exp(-lam*t(i)).*c1*r(i).^(n-2);                              
        else
            bodySource(i) =  -(n^2+lam^2)*exp(-lam*t(i)).*c2*r(i).^(-n-2);
        end
    end
    
end

function uL2 = u_exact_L2_norm(ro,c1,c2,n,lam,halfL)
    uL2 = c1^2*0.5/lam*(1-exp(-lam*2.0*pi))*ro^(2*n+2)/(2*n+2);
    I1 = @(phi) ((halfL*sec(phi)).^(2-2*n)-ro^(2-2*n)).*(exp(-2*lam*phi)+exp(-2*lam*(pi-phi)));
    I2 = @(phi) ((halfL*csc(phi)).^(2-2*n)-ro^(2-2*n)).*exp(-2*lam*phi);
    uL2A = 0.5*c2^2/(1-n)*(quad(I1,0,0.25*pi)+quad(I2,0.25*pi,0.75*pi));
    uL2 = sqrt(uL2+uL2A);
end
