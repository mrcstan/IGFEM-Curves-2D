%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Ahmad Raeisi Najafi on 10/08/2012
%%% Last modified date: 11/02/2013
%%% Copyright 2012 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [UUR, PUR] = solve_matrix_eqn(n_node, KPP,KPF,KFP,KFF,PP,PF,UP, eq_num)

UF =  KFF  \ (PF - KFP * UP);
% Condition = cond(KFF);
% test1 = KFF*UF;
% test2 = (PF - KFP * UP);
% tttt = test1-test2;
% Norm_tttt = norm(tttt, 2);

PP = PP - KPP * UP - KPF * UF; 

% insert the free and prescribed displacements  into Uur
for i = 1:n_node;
    for j = 1:1%j = 1:2
        row = eq_num(i);
        %row = eq_num(i,j);
        if row > 0  % free node
            UUR(i,j) =   UF(row);
            PUR(i,j) =   PF(row);
        elseif row < 0 % prescriibed node
            UUR(i,j) = UP(-row);
            PUR(i,j) = PP(-row);
        end
    end
end

return
