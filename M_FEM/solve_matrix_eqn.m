%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Ahmad Raeisi Najafi on 10/08/2012
%%% Modified by Marcus
%%% Last modified date: 8/21/2014
%%% Copyright 2012 University of Illinois at Urbana-Champaign. All rights reserved
%%% OUTPUT:
%%%     UUR: solution at all nodes and the Lagrange
%%%          multipliers (the last entries after all the dofs)
%%%     PUR: load vector at all nodes and the RHS of the constraint
%%%          equation for Lagrange multipliers
%%%     fGloInd: logical indices of the free nodes and the Lagrange
%%%              multiplier
%%%     pGloInd: logical indices of the prescribed nodes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [UUR, PUR, fGloInd,pGloInd] = solve_matrix_eqn(KPP,KPF,KFP,KFF,PP,PF,UP, eq_num)

UF =  KFF  \ (PF - KFP * UP);
% Condition = cond(KFF);
% test1 = KFF*UF;
% test2 = (PF - KFP * UP);
% tttt = test1-test2;
% Norm_tttt = norm(tttt, 2);

PP = PP - KPP * UP - KPF * UF; 

UUR = zeros(numel(eq_num),1);
PUR = zeros(numel(eq_num),1);
fGloInd = eq_num > 0;
pGloInd = ~fGloInd;
UUR(fGloInd) = UF(eq_num(fGloInd));
UUR(pGloInd) = UP(-eq_num(pGloInd));
PUR(fGloInd) = PF(eq_num(fGloInd));
PUR(pGloInd) = PP(-eq_num(pGloInd));
return
