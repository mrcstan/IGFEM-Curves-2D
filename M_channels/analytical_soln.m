%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Marcus Tan on 1/26/2016
%%% Copyright 2014 University of Illinois at Urbana-Champaign. 
%%% All rights reserved
%%% Analytical solutions for error analysis/convergence study
%%% See "A NURBS-based interface-enriched generalized finite element scheme
%%% for the thermal analysis and design of microvascular composites". CMME,
%%% 283 (2015) 1382-1400 for analytical solutions to a semicircular channel
%%% and a cross network of channels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [u,u_x,u_y] = analytical_soln(x,y)
%% straight channel
% soln can be found in create_channels.m

%% semicircular channel
n = 3;
ro = 0.04;
xc = 0.05;
yc = 0;
% IMPORTANT: Make sure that the mass flow rate and heat capacity from the 
%   channel file, and mesh.material.conductivity are used
conductivity = 1.0;
mcf = 90; % mass flow rate x heat capacity
alpha = mcf/conductivity; % mass flow rate x heat capacity / conductivity
c1 = 1/ro^n;
lam = 2*n/alpha;
c2 = ro^n;
u = u_func(x,y,xc,yc,ro,c1,c2,n,lam);
if nargout > 1
    u_x = ux_func(x,y,xc,yc,ro,c1,c2,n,lam);
    u_y = uy_func(x,y,xc,yc,ro,c1,c2,n,lam);
end

%% cross channel
% soln can be found in create_channels.m or paper
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


