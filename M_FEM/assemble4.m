%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Ahmad Raeisi Najafi on 10/08/2012
%%% Modified by Marcus Tan
%%% Last modified date: 8/7/2014
%%% Copyright 2012 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this function assembles the global stiffness matrix K and load vector P from
% the element stiffness matrices KEL and load vectors PEL
% Only blocks of K, i.e., KPP, KPF, KFP,KFF, where K= [KPP,KPF;KFP,KFF] are
% assembled. Similarly, PP and PF of P=[PP,PF] are assembled.
% The non-zero values are first stored in the form of triple arrays, where the
% first array contains the row indices, the second array contains the column indices and
% the third array contains the values
% For example, for KFF, we have Iff, Jff and valff
% The sparse matrices KPP, KPF, KFP and KFF are then assembled in a single
% call to sparse. For example, KFF = sparse(Iff,Jff,valff,n_dof,n_dof)

function [KPP,KPF,KFP,KFF,PP,PF,UP]=assemble(mesh,eq_num,n_dof,channels,gauss,supg)

% Assemble KPP,KPF,KFP,KFF,PF and UP


UP = assemble_prescribed_node(mesh.node.Dirichlet,eq_num);
%
PP = zeros(mesh.node.Dirichlet.n_pre_temp, 1);
PF = zeros(n_dof, 1);


PF = force_assemble(mesh.elem.elem_node,mesh.node.coords, mesh.elem.Neumann, eq_num, PF);

% preallocate arrays
[estnff,estnfp,estnpp] = estimate_nff_nfp_npp(mesh.elem.n_elem, ...
                                              size(mesh.elem.elem_node,2), ...
                                              mesh.elem.nIGFEMelems, ...
                                              2, ... % ave num of boundary elements sharing one node
                                              2, ... % number of nodes per edge
                                              mesh.node.Dirichlet.n_pre_temp);
fprintf('\nestimated lengths (nff,nfp,npp) = (%i, %i, %i) \n',estnff,estnfp,estnpp)                                           
nff = 0;
Iff = zeros(estnff,1);
Jff = zeros(estnff,1);
valff = zeros(estnff,1);


nfp = 0;
Ifp = zeros(estnfp,1);
Jfp = zeros(estnfp,1);
valfp = zeros(estnfp,1);

npf = 0;
Ipf = zeros(estnfp,1);
Jpf = zeros(estnfp,1);
valpf = zeros(estnfp,1);

npp = 0;
Ipp = zeros(estnpp,1);
Jpp = zeros(estnpp,1);
valpp = zeros(estnpp,1);

% arrays are grown when needed
%{
Iff = [];
Jff = [];
valff = [];
Ifp = [];
Jfp = [];
valfp = [];
Ipf = [];
Jpf = [];
valpf = [];
Ipp = [];
Jpp = [];
valpp = [];
%}
%fprintf('\nBefore building sparse matrix\n')
for i = 1:mesh.elem.n_elem 
    if (mesh.elem.parent(i).type > 1)
        
        [KEL,PEL,errflag] ...
            =compute_parent_element_nurbs(mesh.node.coords,...
                                          mesh.elem.parent(i),...
                                          mesh.elem.heatSource(i),...
                                          mesh.heatSourceFunc,...
                                          mesh.convect,...
                                          channels,...
                                          gauss,...
                                          supg);
       
        % pass the node conectivity of parent element to the current element node conectivity matrix  
        if (any(errflag))
            warning(['Jacobian is non-positive elem num: ',num2str(i),', child(s): ',num2str(find(errflag))])
        end
        curElemNode = [mesh.elem.parent(i).nodes,mesh.elem.parent(i).cstrRows]; 

    else
        %regularTimer = tic;
        [KEL,PEL,errflag] ...
            =compute_regular_element(mesh.node.coords,...
                                     mesh.elem.elem_node(i,:),...
                                     mesh.elem.parent(i),...
                                     mesh.elem.heatSource(i),...
                                     mesh.material(mesh.elem.material(i)).conductivity,...
                                     mesh.heatSourceFunc,...
                                     mesh.convect,...
                                     channels,gauss,supg); 
         % pass the node conectivity of uncut element to the current element node conectivity matrix
         if (errflag)
             warning(['Jacobian is non-positive for elem num: ',num2str(i)]);
         end
         curElemNode = mesh.elem.elem_node(i,:);
         %totRegularTime = totRegularTime + toc(regularTimer);
    end
    
    gloInds = eq_num(curElemNode);
    pLocInds = gloInds < 0;
    fLocInds = gloInds > 0;
    pGloInds = -gloInds(pLocInds);
    fGloInds = gloInds(fLocInds);
    
    % arrays have been preallocated
    % KFF
    [Jel,Iel] = meshgrid(fGloInds,fGloInds);
    nnew = numel(Iel)+nff;
    Iff(nff+1:nnew) = Iel(:);
    Jff(nff+1:nnew) = Jel(:);
    KELsub = KEL(fLocInds,fLocInds); 
    valff(nff+1:nnew) = KELsub(:);
    nff = nnew;

    % KFP
    [Jel,Iel] = meshgrid(pGloInds,fGloInds);
    nnew = numel(Iel)+nfp;
    Ifp(nfp+1:nnew) = Iel(:);
    Jfp(nfp+1:nnew) = Jel(:);
    KELsub = KEL(fLocInds,pLocInds); 
    valfp(nfp+1:nnew) = KELsub(:);
    nfp = nnew;
    
    % KPF    
    [Jel,Iel] = meshgrid(fGloInds,pGloInds);
    nnew = numel(Iel)+npf;
    Ipf(npf+1:nnew) = Iel(:);
    Jpf(npf+1:nnew) = Jel(:);
    KELsub = KEL(pLocInds,fLocInds); 
    valpf(npf+1:nnew) = KELsub(:);
    npf = nnew;
    
    % KPP
    [Jel,Iel] = meshgrid(pGloInds,pGloInds);
    nnew = numel(Iel)+npp;
    Ipp(npp+1:nnew) = Iel(:);
    Jpp(npp+1:nnew) = Jel(:);
    KELsub = KEL(pLocInds,pLocInds); 
    valpp(npp+1:nnew) = KELsub(:);
    npp = nnew;
    
    % the arrays are grown when needed
    %{
    % KFF
    [Jel,Iel] = meshgrid(fGloInds,fGloInds);
    Iff = [Iff;Iel(:)];
    Jff = [Jff;Jel(:)];
    KELsub = KEL(fLocInds,fLocInds); 
    valff = [valff;KELsub(:)];
    
    % KFP
    [Jel,Iel] = meshgrid(pGloInds,fGloInds);
    Ifp = [Ifp;Iel(:)];
    Jfp = [Jfp;Jel(:)];
    KELsub = KEL(fLocInds,pLocInds); 
    valfp = [valfp;KELsub(:)];
    
    % KPF
    [Jel,Iel] = meshgrid(fGloInds,pGloInds);
    Ipf = [Ipf;Iel(:)];
    Jpf = [Jpf;Jel(:)];
    KELsub = KEL(pLocInds,fLocInds); 
    valpf = [valpf;KELsub(:)];
    
    % KPP
    [Jel,Iel] = meshgrid(pGloInds,pGloInds);
    Ipp = [Ipp;Iel(:)];
    Jpp = [Jpp;Jel(:)];
    KELsub = KEL(pLocInds,pLocInds); 
    valpp = [valpp;KELsub(:)];
    %}
    PP(pGloInds) = PP(pGloInds)+PEL(pLocInds);
    
    PF(fGloInds) = PF(fGloInds)+PEL(fLocInds);
    
                                                   
end
%toc

% form each block stiffness matrice with a single call to sparse
fprintf('\nactual lengths (nff,nfp,npp) = (%i, %i, %i) \n',nff,nfp,npp)
KFF = sparse(Iff(1:nff),Jff(1:nff),valff(1:nff),n_dof,n_dof);
KFP = sparse(Ifp(1:nfp),Jfp(1:nfp),valfp(1:nfp),n_dof,mesh.node.Dirichlet.n_pre_temp);
KPF = sparse(Ipf(1:npf),Jpf(1:npf),valpf(1:npf),mesh.node.Dirichlet.n_pre_temp,n_dof);
KPP = sparse(Ipp(1:npp),Jpp(1:npp),valpp(1:npp),mesh.node.Dirichlet.n_pre_temp,mesh.node.Dirichlet.n_pre_temp);

%fprintf('IGFEM timings \n')
%disp(IGFEMtimings)

%fprintf('total IGFEM time = %g \n', totIGFEMtime)
%fprintf('total regular time = %g \n', totRegularTime)
%{
nff = length(Iff)
nfp = length(Ifp)
npf = length(Ipf)
npp = length(Ipp)

KFF = sparse(Iff,Jff,valff,n_dof,n_dof);
KFP = sparse(Ifp,Jfp,valfp,n_dof,mesh.node.Dirichlet.n_pre_temp);
KPF = sparse(Ipf,Jpf,valpf,mesh.node.Dirichlet.n_pre_temp,n_dof);
KPP = sparse(Ipp,Jpp,valpp,mesh.node.Dirichlet.n_pre_temp,mesh.node.Dirichlet.n_pre_temp);
%}
end