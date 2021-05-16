close all
clear
path(path,'../../MatlabUsefulFunctions/export_fig')
ftsz = 36;
nx = 200; % number of grid points in x-direction - 1
ny = 200; % number of grid points in x-direction - 1
xi = 0.0; % left boundary of domain
xf = 0.05; % right boundary of domain
yi = 0.0; % bottom boundary of domain
yf = 0.05; % bottom boundary of domain
x = linspace(xi,xf,nx+1);
y = linspace(yi,yf,ny+1);

[XX,YY] = meshgrid(x,y);

% For validation
%{
% [xi,xf,yi,yf] = [0,0.15,0,0.2]
func = @(x,y) 500*(0.2489 + 21.87*y + 23.42*x - 356.8*y.^2 - 16.96*x.*y ...
            -485.7*x.^2 + 2422*y.^3 + 276.7*y.^2.*x + 170.4*y.*x.^2 ...
            +4317*x.^3 - 5849*y.^4 - 972.2*y.^3.*x - 490*x.^2.*y.^2 ...
            -791*y.*x.^3 - 14060*x.^4);
val = func(XX,YY);
%}

% Single localized source
%{
% [xi,xf,yi,yf] = [0,0.15,0,0.2]
xo = 0.075;
yo = 0.1;
ro = 0.04;
qo = 500*0.15*0.2*(15/(16*ro))^2;
func = @(x,y) ((x>=(xo-ro)) & (x<=(xo+ro)) & (y>=(yo-ro)) & (y<=(yo+ro)))...
               .*qo.*(1-((x-xo)/ro).^2).^2.*(1-((y-yo)/ro).^2).^2 ...
               + ((x<(xo-ro)) | (x>(xo+ro)) | (y<(yo-ro)) | (y<(yo+ro)))*0; 
val = func(XX,YY)/1000.0;
%}

% Double localized sources
%{
% [xi,xf,yi,yf] = [0,0.15,0,0.2]
xo = 0.04;
yo = 0.04;
ro = 0.015;
qo = 250*0.15*0.2*(15/(16*ro))^2;
x1 = 0.11;
y1 = 0.16;
r1 = 0.015;
q1 = 250*0.15*0.2*(15/(16*ro))^2;
func = @(x,y) ((x>=(xo-ro)) & (x<=(xo+ro)) & (y>=(yo-ro)) & (y<=(yo+ro)))...
               .*qo.*(1-((x-xo)/ro).^2).^2.*(1-((y-yo)/ro).^2).^2 ...
               + ((x>=(x1-r1)) & (x<=(x1+r1)) & (y>=(y1-r1)) & (y<=(y1+r1)))...
               .*q1.*(1-((x-x1)/r1).^2).^2.*(1-((y-y1)/r1).^2).^2; 
val = func(XX,YY)/1000.0;
%}

% For grid-like network, redundancy/damage-resilience project
% in meters
%
func = @(x,y) 6000*(0.188625 + 127*x + 99.54*y ...
                  -7.583E3*x.^2 - 4.038E3*x.*y - 5.169E3*y.^2 ...
                  +1.919E5*x.^3 + 2.205E5*x.^2.*y + 6.490E4*y.^2.*x ...
                  +1.228E5*y.^3 - 1.700E6*x.^4 - 5.539E6*x.^3.*y ...
                  -1.921E6*x.^2.*y.^2 + 1.263E6*x.*y.^3 - 2.156E6*y.^4 ...
                  -2.348E6*x.^5 + 5.443E7*x.^4.*y - 8.962E6*x.^3.*y.^2 ...
                  +4.017E7*x.^2.*y.^3 - 4.462E7*x.*y.^4 + 1.982E7*y.^5);
%
% in mm
% [xi,xf,yi,yf] = [0,50,0,50]
%{
func = @(x,y) 6000*(0.188625 + 0.1270*x + 9.954E-02*y ...
                   -7.583E-03*x.^2 - 4.038E-03*x.*y - 5.169E-03*y.^2 ...
                   +1.919E-04*x.^3 + 2.205E-04*x.^2.*y + 6.490E-05*y.^2.*x ...
                   +1.228E-04*y.^3 - 1.700E-06*x.^4 - 5.539E-06*x.^3.*y ...
                   -1.921E-06*x.^2.*y.^2 + 1.263E-06*x*y.^3 - 2.156E-06*y.^4 ...
                   -2.348E-09*x.^5 + 5.443E-08*x.^4.*y - 8.962E-09*x.^3.*y.^2 ...
                   +4.017E-08*x.^2.*y.^3 - 4.462E-08*x.*y.^4 + 1.982E-08*y.^5);              
%}
val = func(XX,YY);

%surf(XX,YY,val,'edgecolor','none')
%shading interp
contourf(XX,YY,val)
axis image
h = colorbar;
set(h,'fontsize',ftsz)
ylabel(h,'$f(x,y)\,\mathrm{(kW/m^2)}$ ','Interpreter','LaTex')
xlabel('$x$ (m)','Interpreter','LaTex','fontsize',ftsz)
ylabel('$y$ (m)','Interpreter','Latex','fontsize',ftsz)
%xlim([0,0.151])
set(gca,'fontsize',ftsz)
%set(gca,'units','normalized','position',[0.15,0.15,0.7,0.7],'xtick',[0,0.05,0.1,0.15])
