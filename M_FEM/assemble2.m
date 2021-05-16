%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Ahmad Raeisi Najafi on 10/08/2012
%%% Modified by Marcus Tan
%%% Last modified date: 09/30/2013
%%% Copyright 2012 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [KPP,KPF,KFP,KFF,PP,PF,UP]=assemble(mesh,eq_num,n_dof,itrface,gauss,supg)

% Assemble KPP,KPF,KFP,KFF,PF and UP


UP = assemble_prescribed_node(mesh.node.Dirichlet,eq_num);
%
PP = zeros(mesh.node.Dirichlet.n_pre_temp, 1);
PF = zeros(n_dof, 1);
%
%{
PP = sparse(mesh.allNode.Dirichlet.n_pre_temp, 1);
PF = sparse(n_dof, 1);
%}

PF = force_assemble(mesh.elem.elem_node,mesh.node.coords, mesh.elem.Neumann, eq_num, PF);

%{
KPP = zeros(mesh.allNode.Dirichlet.n_pre_temp, mesh.allNode.Dirichlet.n_pre_temp);
KPF = zeros(mesh.allNode.Dirichlet.n_pre_temp, n_dof);
KFP = zeros(n_dof, mesh.allNode.Dirichlet.n_pre_temp);
KFF = zeros(n_dof, n_dof);                   
%}                
%
KPP = sparse(double(mesh.node.Dirichlet.n_pre_temp), double(mesh.node.Dirichlet.n_pre_temp));
KPF = sparse(double(mesh.node.Dirichlet.n_pre_temp), double(n_dof));
KFP = sparse(double(n_dof), double(mesh.node.Dirichlet.n_pre_temp));
KFF = sparse(double(n_dof), double(n_dof));
%{
elemSet = [];
gauss2 = gauss;
gauss2.qua1Dt = gauss_points_and_weights(false,2,1); % normal to channel
gauss2.qua1Dn = gauss_points_and_weights(false,2,1); % normal to channel
%}

disp('assembly of stiffness matrix')
tic 
for elem_num = 1:mesh.elem.n_elem 

    if (mesh.elem.parent(elem_num).isParent)
        %{
        if(any(elemSet == elem_num))
            disp('reduced integration for elem:')
            disp(elem_num)
             [KEL,PEL]=compute_parent_element_nurbs(mesh.node.coords,...
                mesh.elem.parent(elem_num),...
                itrface,mesh.heatSource,gauss2);
        else
        %}   
        [KEL,PEL,errflag] ...
            =compute_parent_element_nurbs(mesh.node.coords,...
                                          mesh.elem.parent(elem_num),...
                                          mesh.elem.heatSource(elem_num),...
                                          mesh.heatSourceFunc,...
                                          mesh.convect,...
                                          itrface,...
                                          gauss,...
                                          supg);
        %end
        % pass the node conectivity of parent element to the current element node conectivity matrix  
        if (any(errflag))
            warning(['Jacobian is non-positive elem num: ',num2str(elem_num),', child(s): ',num2str(find(errflag))])
        end
        current_elem_node = mesh.elem.parent(elem_num).node; 
    else
        [KEL,PEL,errflag] ...
            =compute_regular_element(mesh.node.coords,...
                                     mesh.elem.elem_node(elem_num,:),...
                                     mesh.elem.parent(elem_num),...
                                     mesh.elem.heatSource(elem_num),...
                                     mesh.material(mesh.elem.material(elem_num)).conductivity,...
                                     mesh.heatSourceFunc,...
                                     mesh.convect,...
                                     itrface,gauss,supg); 
         % pass the node conectivity of uncut element to the current element node conectivity matrix
         if (errflag)
             warning(['Jacobian is non-positive for elem num: ',num2str(elem_num)]);
         end
         current_elem_node = mesh.elem.elem_node(elem_num,:);
    end
                          
    [KPP,KPF,KFP,KFF,PP,PF] = element_stiffness_assemble2 ...
        (current_elem_node,eq_num,KEL, PEL, ... 
         KPP,KPF,KFP,KFF,PP,PF);
end
toc
return