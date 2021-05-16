%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Ahmad Raeisi Najafi on 10/08/2012
%%% Last modified date: 01/03/2013
%%% Copyright 2012 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function maps the global node number to the equation number,
% with the free nodes mapped to positive equation number and prescribed 
% nodes mapped to negative equation number
function [n_dof, eq_num] = initialize(n_node, n_pre_temp, temp_node)



% Initialize the equation number vector
eq_num  = zeros(n_node, 1);

% assign "negative" equation numbers for nodes with prescribed displacements
for i = 1:n_pre_temp
    node = temp_node(i);
    %idof = disp_node(i, 2);
    eq_num(node) = -i;
end

% assign equation numbers for free nodes
n_dof = 0;
for j = 1:n_node
    %for m = 1:2
        if eq_num(j) == 0
            n_dof = n_dof + 1;
            eq_num(j) = n_dof;
        end
    %end
end
% assign equation numbers for element corners

return; % end of function