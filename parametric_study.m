 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Marcus Tan
%%% Modified by Marcus Tan on 1/17/2015
%%% Copyright 2012 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%function Main
%%%%%%%%%%%%%%%%%%%%%%%%%
%  IGFEM Code  
%%%%%%%%%%%%%%%%%%%%%%%%%
%
clear all
%
close all
format long e;

% paths
path(path, './M_geom_toolbox')
path(path, './M_preFEM')
path(path, './M_channels')
path(path, './M_FEM')
path(path, './mx_FEM')
path(path, './M_postprocessing')
path(path, './M_error_analysis')
path(path, '../NURBS/nurbs_toolbox')
path(path, './mesh_conforming_abaqus')
path(path, './mesh')
path(path, '../SISL') % for windows
path(path, './ChannelFiles')
path(path, '../Opt-IGFEM-Curves-2D/M_optimization')
path(path, '../MatlabUsefulFunctions/export_fig')

totTimer = tic;
%% MESH AND USER INPUT
% 2 choices: I) supply an Abaqus mesh file
%            II) give information for generating the mesh in this file
% In both cases, the BC's must be specify in this file as this code is
% capable of mesh refinement
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% II) Information for generating mesh
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% domain boundary
% 
% conforming mesh for parallel channels with nx=15, ny=18
%{
mesh.boundary.xi = 0.0;
mesh.boundary.xf = 0.15;
mesh.boundary.yi = -0.00125;
mesh.boundary.yf = 0.20125;
%}
%
mesh.boundary.xi = 0.0;
mesh.boundary.xf = 0.15;
mesh.boundary.yi = 0.0;
mesh.boundary.yf = 0.2;
%
 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Other user inputs that must be specified regardless of choice I or II
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Interface, analytical body source (if available) and analytical 
% solution (if available)
svfile = 'parallel2_movable_mid_channel_xo_p04_yo_p1_x1_p11_y1_p1_ror1_p015_QoQ1_250_npts_21x21_30x40_mesh';
channelFile = 'parallel2_movable_mid_channel.channel';
[channels,mesh.heatSourceFunc,designParams] = preprocess_channels(channelFile);
% boundary conditions
% 1 : x = xi, 2 : x = xf, 3 : y = yi, 4 : x = yf
mesh.BCs.boundaries = []; % totally insulated boundaries
%mesh.BCs.boundaries = [3,4];
%mesh.BCs.types = [2,1];
%mesh.BCs.values_or_funcs = {30,27};
%mesh.BCs.boundaries = [3,4]; 
%mesh.BCs.types = [2,1];
%mesh.BCs.values_or_funcs = {2000,20};
%mesh.BCs.values_or_funcs = {2000,@(x,y) 20*x/0.1};
%mesh.BCs.values_or_funcs = {1,21.5};
%mesh.BCs.boundaries = [1,2,3,4];
%mesh.BCs.types = [1,1,1,1];
%mesh.BCs.values_or_funcs = {u,u,u,u};   


% specify convective heat coefficient 
%mesh.convect.coef = 8.5+5.71; old value
%mesh.convect.coef = 8+5.6;
%mesh.convect.Tref = 21;
mesh.convect.coef = 0.0; 
mesh.convect.Tref = 0.0;
  

% for nodes_curves_intersect.m and move_nodes.m 
tol.node = 1e-6; % tolerance for checking if original nodes coincide with channels
moveNode.distFrac = 0.05; % distance to move node when it coincides with channel
                         % or when a branching point or kink coincides
                         % with an element edge, as a fraction of the
                         % minimum edge length
moveNode.maxAttemps = 5;
tol.boundary = 1e-13; % also for set_boundary_conditions.m

% tolerances for edges_curves_intersect.m
tol.nurbsParam = 1.e-8;
tol.epsco = 1e-15; % computational tolerance for curve_edge_intersect.mexw64
tol.epsge = 1e-6; % geometrical tolerance for curve_edge_intersect.mexw64
tol.cosAngleTol = cos(20*pi/180.0); % cosine of the min angle of the resulting element after edge flipping
tol.halfLineWidthFrac = 1e-4; % half-width fraction of interface
tol.vert = tan(89.9*pi/180); % tan of max positive angle wrt horizontal axis beyond which a line is considered vertical
opt.maxRefineLevel = 12;
opt.refineJuncElem = false; % force refinement of element with branching or kinks
tol.intersectEdges = 1e-13; % for intersect_edges in single_edge_curves_intersect.m


% tolerenace for approximating slender child element with triangle or quadrilateral
% slenderTol = []; % if approximation not desired
slenderTol.minAngle = 5; % min angle of a triangular child element below which NURBS is approximated by line segments
slenderTol.maxAspectRatio = inf;  % max aspect ratio of a quadrialteral child element below which NURBS is approximated by line segments


% integration schemes
triNpt1D = 4; % for line integration in regular FEM,poly IGFEM and NURBS IGFEM
triNpt2D = 16; % for element integration in regular FEM and poly IGFEM
quaNpt1Dt = 4; % number of gauss points per knot span along channel for NURBS IGFEM 
quaNpt1Dn = 4; % number of gauss points per knot span normal to channel for NURBS IGFEM
quadNpt1D = 4; % number of gauss points in one direction for quad child element in polynomial IGFEM

    
% perform FEM analysis
performFEM = true;
polyIGFEM = true;
supg = true; % apply SUPG. Note: do not apply SUPG for constant heat flux model

nParamPts = [21,21]; % number of equally space points for each design parameter within its bound


%% Calculate solution for each user defined point in design parameter space 
totParamPts = prod(nParamPts);
if (designParams.nParams == 1)
    designParamPts = linspace(designParams.bounds(1,1), ...
                              designParams.bounds(1,2), ...
                              nParamPts)';
    objVals = nan(totParamPts); 
    objVals2 = nan(totParamPts); 
elseif (designParams.nParams == 2)
    [designParamPt1,designParamPt2] ...
        = meshgrid(linspace(designParams.bounds(1,1), ...
                            designParams.bounds(1,2), ...
                            nParamPts(1)), ...
                   linspace(designParams.bounds(2,1), ...
                              designParams.bounds(2,2), ...
                              nParamPts(2)));
    designParamPts = reshape(cat(2,designParamPt1,designParamPt2),[],2);
    objVals = zeros(nParamPts(1),nParamPts(2));
    objVals2 = zeros(nParamPts(1),nParamPts(2));
else
    error('more than two design parameters not considered')
end


figure(7)
set(gcf,'units','normalized','outerposition',[0 0 1 1])
contourf(designParamPt1,designParamPt2,objVals,12)
xlabel('d_1 (m)','fontsize',30)
ylabel('d_2 (m)','fontsize',30)
set(gca,'fontsize',30)
hcolor = colorbar;
set(hcolor,'fontsize',30)
ylabel(hcolor,'T_8 (^{\circ}C)','rot',0)
axis image
    
    
for sim = 1:totParamPts    
    fprintf('\n---------------------------------------------------------\n')
    fprintf('Simulation %i of %i',sim, totParamPts)
    fprintf('\n---------------------------------------------------------\n')
    % element connectivity and node coordinates
    [mesh.elem.elem_node,mesh.node.coords] ...
        = generate_uniform_mesh([30,40],... % number of elements in x- and y-directions
                                [mesh.boundary.xi,...
                                 mesh.boundary.xf,...
                                 mesh.boundary.yi,...
                                 mesh.boundary.yf],...
                                 2); % 1: diagonal along NW direction
                                     % 2: diagonal along NE direction


    % element material and material conductivity
    %mesh.material.conductivity = 0.2*0.003; % epoxy
    %mesh.material.conductivity = 0.00288; % epoxy
    %mesh.material.conductivity = 2.9*0.003; % old composite value
    mesh.material.conductivity = 2.7*0.003; % composite value Dec 16, 2014
    %mesh.material.conductivity = 0.6;
    %mesh.material.conductivity = 1.0;
    mesh.elem.material = int32(ones(size(mesh.elem.elem_node,1),1));

    % distributed heat source 
    mesh.elem.heatSource = 0*ones(size(mesh.elem.elem_node,1),1);
    
    %stats = mesh_statistics(mesh.elem.elem_node,mesh.node.coords);
    % label edges for refinement later
    mesh.elem.elem_node = label(mesh.node.coords,mesh.elem.elem_node); 

    % the edge_node information is generated by generate_conforming_mesh
    fprintf('\ngenerate edge_node\n')
    tic
    mesh.edge.edge_node = find_edge_node(mesh.elem.elem_node);
    toc

    mesh.edge.length = find_edge_length(mesh.edge.edge_node,mesh.node.coords);

    mesh.edge.minLength = min(mesh.edge.length);
    tol.halfLineWidth = mesh.edge.minLength*tol.halfLineWidthFrac;
    moveNode.dist = mesh.edge.minLength*moveNode.distFrac;
    channels= update_channels(designParamPts(sim,:), ...
                              designParams, ...
                              [], ...
                              channels, ...
                              'replace',...
                              false);
    all_figs = findobj(0, 'type', 'figure');
    close(setdiff(all_figs,7)); % close all figures except the history figure                      
    %plot_mesh_curve(mesh.node.coords, ...
    %                mesh.elem.elem_node, ...
    %                channels,false,false)
    % find intersection points and enrichment nodes
    fprintf('\nfinding intersection points \n')
    tic    
    calcItrsectVel = false;
    [mesh.edge,mesh.node,mesh.elem]...
         =edges_curves_intersect(mesh.edge,mesh.node,...
                                 mesh.elem, ...
                                 channels, ...
                                 tol, ...
                                 opt, ...
                                 calcItrsectVel);


    toc
    

    %% Set boundary conditions
    fprintf('\nsetting boundary conditions \n')
    tic
    [mesh.elem, mesh.node] = set_boundary_conditions(mesh.BCs,...
                                                     mesh.elem,...
                                                     mesh.node,...
                                                     mesh.boundary,...
                                                     tol.boundary);
    toc

    %% enrichment functions and integration subdomains
    % creat parent elements
    fprintf('\nconstructing parent elements\n')
    tic
    mesh.elem.parent = element_intersections(mesh.elem.elem_node, ...
                                             mesh.elem.dualedge, ...
                                             mesh.edge.edge_node, ...
                                             mesh.edge.itrsect, ...
                                             calcItrsectVel);
    %
    [mesh.elem.parent, ...
     mesh.elem.cstrElems, ...
     mesh.elem.nIGFEMelems] ...
                     = parent_elements_nurbs(mesh.elem.parent,...
                                             mesh.elem.elem_node,...
                                             mesh.elem.elem_edge,...
                                             mesh.elem.material,...
                                             mesh.material,...
                                             mesh.elem.junc,...
                                             mesh.elem.kinks,...
                                             mesh.node.coords,...
                                             mesh.node.nOriginalNode,...
                                             mesh.node.constraint,...
                                             channels,...
                                             tol.nurbsParam,...
                                             polyIGFEM,...
                                             calcItrsectVel);


    toc
    elem2plot = [];

    if(mesh.node.n_node==mesh.node.nOriginalNode)
        fprintf('\nconforming mesh\n')
    else
        fprintf('\nnon-conforming mesh\n')
        fprintf('constructing child elements\n')
        tic

        [mesh.elem.parent,mesh.node]...edge, Dirichlet, nodeCoords, boundary
            =child_elements_nurbs(mesh.elem.parent,... 
                                  mesh.node,...
                                  slenderTol,...
                                  polyIGFEM);
        toc                      
        % the global equation number of additional equations for 
        % Lagrange multiplier method                     
        for i = 1:numel(mesh.elem.cstrElems)
            mesh.elem.parent(mesh.elem.cstrElems(i)).cstrRows...
                    =  mesh.elem.parent(mesh.elem.cstrElems(i)).cstrRows ...
                     + size(mesh.node.coords,1); 
        end

    end


    %%  Finite Element Code  
    %%%%%%%%%%%%%%%%%%%%%%%%%  
    % gauss quadrature schemes
    %

    fprintf('\nperforming FEM\n')

    % eq_num:  Equation number assigned to each node
    % n_dof: The number of degree of freedom in the model
    [n_dof, eq_num] = initialize(size(mesh.node.coords,1) ...
                                 +numel(mesh.node.constraint.temp_node), ...
                                 mesh.node.Dirichlet.n_pre_temp, ...
                                 mesh.node.Dirichlet.temp_node);

    % Assemble the stiffness matrix
    gauss.line = gauss_points_and_weights(true,triNpt1D,1,'combined');
    gauss.elem = gauss_points_and_weights(true,triNpt2D ,2,'combined');
    timer = tic;
    if (polyIGFEM)
        fprintf('\nassembly of polynomial IGFEM stiffness matrix\n')
        gauss.quadElem = gauss_points_and_weights(false,quadNpt1D,2,'combined');
        %
        UP =  assemble_prescribed_node(mesh.node.Dirichlet,eq_num);
        [KFF,KFP,KPF,KPP,PF,PP] ...
                    = mx_assemble_sparse(mesh.node.coords',...
                                         mesh.elem.elem_node',...
                                         mesh.elem.heatSource,...
                                         mesh.convect,...
                                         eq_num,...
                                         gauss,...
                                         mesh.elem.parent,...
                                         channels,...
                                         mesh.node.Dirichlet,...
                                         mesh.elem.Neumann,...
                                         supg);
        
    else
        gauss.qua1Dt = gauss_points_and_weights(false,quaNpt1Dt,1,'separate'); % along channel
        gauss.qua1Dn = gauss_points_and_weights(false,quaNpt1Dn,1,'separate'); % normal to channel
        fprintf('assembly \n')
        [KPP,KPF,KFP,KFF,PP,PF,UP] = ...
            assemble (mesh, eq_num, n_dof,channels,gauss,supg);
    end
    toc(timer)
    
    
    
    fprintf('\nsolving equation\n')
    tic
    [UUR, PUR] = solve_matrix_eqn(KPP,KPF,KFP,KFF,PP,PF,UP, eq_num);  
    toc            
    
    % output paraview file
    fprintf('\nupdate enrichment node values\n')
    tic
    UUR2 = update_enrichment_node_value(UUR,mesh.node,mesh.elem,mesh.edge);
    toc
    
    fprintf('\n total simulation time %g \n',toc(totTimer))
    
    %fprintf('\n average temp = %g \n',average_temp(mesh.elem,mesh.node.coords,UUR2));
    Tmax = max(UUR2);
    fprintf('\n max temp = %g \n',Tmax);
    normp = 8;
    objVals(sim) = field_p_norm(mesh.elem,mesh.node.coords,UUR,...
                                normp,Tmax,gauss);
    normp = 16;
    objVals2(sim) = field_p_norm(mesh.elem,mesh.node.coords,UUR,...
                                normp,Tmax,gauss);                        
    figure(7)
    contourf(designParamPt1,designParamPt2,objVals,12)
    xlabel('d_1 (m)','fontsize',30)
    ylabel('d_2 (m)','fontsize',30)
    set(gca,'fontsize',30)
    hcolor = colorbar;
    set(hcolor,'fontsize',30)
    ylabel(hcolor,'T_8 (^{\circ}C)','rot',0)
    axis image
    
end
objVals = objVals/(0.15*0.2)^(1/8);
objVals2 = objVals2/(0.15*0.2)^(1/16);
save(svfile)
figure(7)
contourf(designParamPt1,designParamPt2,objVals)
xlabel('d_1 (m)','fontsize',30)
ylabel('d_2 (m)','fontsize',30)
set(gca,'fontsize',30)
hcolor = colorbar;
set(hcolor,'fontsize',30)
ylabel(hcolor,'T_8 (^{\circ}C)','rot',0)
axis image

figure(8)
surf(designParamPt1,designParamPt2,objVals)
shading interp
xlabel('d_1 (m)','fontsize',30)
ylabel('d_2 (m)','fontsize',30)
zlabel('T_8 (^{\circ}C)','fontsize',30)
view(30,30)
set(gca,'fontsize',30)

figure(9)
contourf(designParamPt1,designParamPt2,objVals2,12)
xlabel('d_1 (m)','fontsize',30)
ylabel('d_2 (m)','fontsize',30)
set(gca,'fontsize',30)
hcolor = colorbar;
set(hcolor,'fontsize',30)
ylabel(hcolor,'T_{16} (^{\circ}C)','rot',0)
axis image

figure(10)
surf(designParamPt1,designParamPt2,objVals2)
shading interp
xlabel('d_1 (m)','fontsize',30)
ylabel('d_2 (m)','fontsize',30)
zlabel('T_{16} (^{\circ}C)','fontsize',30)
view(30,30)
set(gca,'fontsize',30)

fclose ('all');
%