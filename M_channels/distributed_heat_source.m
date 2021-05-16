%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Marcus Tan on 1/26/2016
%%% Copyright 2014 University of Illinois at Urbana-Champaign. 
%%% All rights reserved
%%% Analytical solutions for error analysis/convergence study
%%% See "A NURBS-based interface-enriched generalized finite element scheme
%%% for the thermal analysis and design of microvascular composites". CMME,
%%% 283 (2015) 1382-1400 for analytical solutions to a semicircular channel
%%% and a cross network of channels
%%% NOTE: The distributed_heat_source defined here is only applicable for 
%%%       NURBS-based IGFEM.
%%%       To define one for polynomial IGFEM, one must use the file
%%%       body_source_functions.cpp in the mx_FEM directory
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function val = distributed_heat_source(x,y)
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
lam = 2*n/alpha;
c1 = 1/ro^n;
c2 = ro^n;
val = semicircular_channel_source(x,y,xc,yc,ro,c1,c2,n,lam)/conductivity;
end


function vals = semicircular_channel_source(x,y,xc,yc,ro,c1,c2,n,lam)
    [t,r] = cart2pol(x-xc,y-yc);
    len = length(t);
    vals = zeros(len,1);
    for i = 1:len
        if (r(i) < ro)
            vals(i) = -(n^2+lam^2)*c1*exp(-lam*t(i)).*r(i).^(n-2);                              
        else
            vals(i) =  -(n^2+lam^2)*c2*exp(-lam*t(i)).*r(i).^(-n-2);
        end
    end
    
end