%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Ahmad Raeisi Najafi on 10/08/2012
%%% Modified by Marcus Tan
%%% Last modified date: 06/13/2013
%%% Copyright 2012 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This function assembles the element stiffness matrices into global
% matricies, partitioned according to the equation numbering.
% INPUT:
% elem_node: the current element node conectivity
function [Kpp,Kpf,Kfp,Kff,Pp,Pf] = element_stiffness_assemble2 ...
        (elem_node,eq_num,Kel,Pel,Kpp,Kpf,Kfp,Kff,Pp,Pf)


nElemNode=size(elem_node,2);


for i = 1:nElemNode 
    node_i = elem_node(i);
    row = eq_num(node_i); % row in global system
    for k = 1:nElemNode  
        node_k = elem_node(k);
        col = eq_num(node_k); % column in global system
        if row < 0 && col < 0
            Kpp(-row,-col) = Kpp(-row,-col) + Kel(i,k);
        elseif row < 0 && col > 0
            Kpf(-row,col) = Kpf(-row,col) + Kel(i,k);
        elseif row > 0 && col < 0
            Kfp(row,-col) = Kfp(row,-col) + Kel(i,k);
        elseif row > 0 && col > 0
            Kff(row,col) = Kff(row,col) + Kel(i,k);
        end
    end
    if row < 0
        Pp(-row) = Pp(-row) + Pel(i);
    elseif row > 0 %&& col < 0
        Pf(row) = Pf(row) + Pel(i);
    end           
end
end % end of function