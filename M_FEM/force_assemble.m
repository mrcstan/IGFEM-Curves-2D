 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Ahmad Raeisi Najafi on 10/08/2012
%%% Modified by Marcus Tan on 9/30/2013      
%%% Last modified date: 8/5/2013  
%%% Copyright 2012 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function PF = force_assemble(elem_node,coords, Neumann, eq_num, PF)

% assemble PF

for i = 1: Neumann.n_heatFlux
    %
    element = Neumann.heatFlux_elem(i);
    switch Neumann.heatFlux_surface(i)
        case 1 
            node1 = elem_node(element, 1);
            node2 = elem_node(element, 2);

        case 2
            node1 = elem_node(element, 2);
            node2 = elem_node(element, 3);

        case 3
            node1 = elem_node(element, 3);
            node2 = elem_node(element, 1);

        otherwise
            error('Unexpected Edge.');
    end
    %
    %node1 = Neumann.heatFlux_endNode(1,i);
    %node2 = Neumann.heatFlux_endNode(2,i);
    d = sqrt((coords(node1,1) - coords(node2,1))^2 + ...
              (coords(node1,2) - coords(node2,2))^2); % The length of edge
    f1 = 0.5 * d * Neumann.heatFlux_value(i);
    f2 = 0.5 * d * Neumann.heatFlux_value(i);
    
    row1 = eq_num(node1);
    row2 = eq_num(node2);
  
    
    if row1 > 0
        PF(row1) = PF(row1) + f1;
    end
    
    if row2 > 0
        PF(row2) = PF(row2) + f2;
    end
end

return;