%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Marcus Tan on 1/26/2016
%%% Copyright 2016 University of Illinois at Urbana-Champaign. 
%%% All rights reserved
%%% Analytical solutions for error analysis/convergence study
%%% See "A NURBS-based interface-enriched generalized finite element scheme
%%% for the thermal analysis and design of microvascular composites". CMME,
%%% 283 (2015) 1382-1400 for analytical solutions to a semicircular channel
%%% and a cross network of channels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function val = soln_L2_norm()
n = 3;
ro = 0.04;
% IMPORTANT: Make sure that this is the same as that of the channel file! 
mcf = 1.747; % mass flow rate x heat capacity
c1 = 1/ro^n;
lam = 2*n/mcf;
c2 = c1*ro^(2*n);
%% Semicircular channel 
val = u_exact_L2_norm(ro,c1,c2,n,lam,0.05);
end


function uL2 = u_exact_L2_norm(ro,c1,c2,n,lam,halfL)
    uL2 = c1^2*0.5/lam*(1-exp(-lam*2.0*pi))*ro^(2*n+2)/(2*n+2);
    I1 = @(phi) ((halfL*sec(phi)).^(2-2*n)-ro^(2-2*n)).*(exp(-2*lam*phi)+exp(-2*lam*(pi-phi)));
    I2 = @(phi) ((halfL*csc(phi)).^(2-2*n)-ro^(2-2*n)).*exp(-2*lam*phi);
    uL2A = 0.5*c2^2/(1-n)*(integral(I1,0,0.25*pi)+integral(I2,0.25*pi,0.75*pi));
    uL2 = sqrt(uL2+uL2A);
end