%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Marcus Tan on 3/17/2014
%%% Last modified date: 3/17/2014
%%% Copyright 2014 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% h can only be evaluated one point at a time
function h = conductance(x,kapf,mcf)
factor = 0.5*pi*kapf/mcf;
betas = [25.6796, 83.8618, 174.167, 296.536,450.947,637.387,855.850];
R1 = [-0.492517,0.395508,-0.345872,0.314047,-0.291252,0.273808,-0.259852];
C = [0.403483,-0.175111,0.105594,-0.0732804,0.0550357,-0.043483,0.035597];
h = 2*pi*kapf./(4*x*factor + 11.0/24.0 + sum(C.*R1.*exp(-betas*factor*x)));

end