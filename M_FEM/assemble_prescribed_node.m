%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Ahmad Raeisi Najafi on 10/08/2012
%%% Modified by Marcus Tan on 7/22/2013
%%% Last modified date: 9/7/2013
%%% Copyright 2012 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This function defines the prescribed nodal values. 
% INPUT:
%   dirichlet.n_pre_temp
%   dirichlet.temp_node
%   dirichlet.temp_value
% OUTPUT:
%   UP: the presribed nodal values. Note that the nodal value at a
%   presribed enrichment node is the difference between presribed value and 
%   the interpolated value based on the values at the two adjacent nodes on
%   the original element. 
function UP = assemble_prescribed_node(dirichlet,eq_num)

UP =zeros(dirichlet.n_pre_temp, 1); % Initilize the matrix of UP to zeros
% assemble UP
for i = 1: dirichlet.n_pre_temp
    node = dirichlet.temp_node(i);
    %dof = disp_node(i,2); 
    %row = -eq_num(node, dof);
    row = -eq_num(node);
    UP(row,1) = dirichlet.temp_value(i);
end

end