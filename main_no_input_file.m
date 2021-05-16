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
clear
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
path(path, '../SISL/mx_SISL') % for windows
path(path, './ChannelFiles')
path(path, '../Opt-IGFEM-2D/M_optimization')
path(path, '../MatlabUsefulFunctions/export_fig')
%path(path, '../../SISL_linux') % for LINUX

totTimer = tic;
%% MESH AND USER INPUT
% 2 choices: I) supply an Abaqus mesh file
%            II) give information for generating the mesh in this file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% I) Abaqus mesh file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% straight/cross channel
%meshfile = 'square_p1_10x10_patran_temp_all.inp';
%meshfile = 'square_p1_80x80C_abq_temp_all.inp';
%meshfile = 'square_p1_10x10_abq_temp_top_heatflux_bottom.inp';

% semircircular channel
%meshfile = 'p1xp05_unstructured4_abq_temp_all.inp';
%meshfile = 'p1xp05_40x20_patran_temp_all.inp';
%meshfile = 'conform0_temp_all_semicircle_ro_p02.inp';

% single wavy/network channel
%meshfile = 'p1xp05_10x5_abq_temp_top_heatflux_bottom.inp';
%meshfile = 'wavy_network_conform2b_mesh.inp';

% 6-channels - validation with FLUENT
% meshfile = 'p2xp15_80x60_mesh_insulated.inp';
% path(path, './SixChannels')
% curved_channel_with_circular_arcs - validation with FLUENT
%meshfile = 'p1xp05_40x20_abq_insulated.inp';
%path(path, './CurvedChannelWithCircularArcs');
% serpentine channel
% meshfile = 'square_p1_80x80_abq_insulated.inp';
%
%meshfile = 'p1xp05_4x2_patran_temp_all.inp';
%mesh = read_input(meshfile);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% II) Information for generating mesh
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% domain boundary
mesh.boundary.xi = 0.0;
mesh.boundary.xf = 0.15;
mesh.boundary.yi = 0.0;
mesh.boundary.yf = 0.2;

meshSizes = [15,20];
% element connectivity and node coordinates
[mesh.elem.elem_node,mesh.node.coords] ...
    = generate_uniform_mesh(meshSizes,... % number of elements in x- and y-directions
                             [mesh.boundary.xi,...
                             mesh.boundary.xf,...
                             mesh.boundary.yi,...
                             mesh.boundary.yf],...
                             2); % 1: diagonal along NW direction
                                 % 2: diagonal along NE direction
                               
% element material and material conductivity
%mesh.material.conductivity = 0.2*0.003; % epoxy
%mesh.material.conductivity = 0.00288; % epoxy
% with convection: 3D conductivity = 2.04
% without convection: 3D conductivity = 2.7
% PDMS conductivity: 0.3
mesh.material.conductivity = 2.7*0.003; % composite value Dec 16, 2014
%mesh.material.conductivity = 1.0;
mesh.elem.material = int32(ones(size(mesh.elem.elem_node,1),1));

% distributed heat source 
mesh.elem.heatSource = 500.0*ones(size(mesh.elem.elem_node,1),1); 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Other user inputs that must be specified regardless of choice I or II
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Interface, analytical body source (if available) and analytical 
% solution (if available)
channelFile = 'parallel2_start.channel';
%channelFile = 'DK_DP_check.channel';
options.suppressDesignParams = false;
[channels,designParams] = preprocess_channels(channelFile,options);

% calculate volume fraction
%domainVol = (mesh.boundary.xf - mesh.boundary.xi)*(mesh.boundary.yf - mesh.boundary.yi)*0.003; 
%fprintf('volume fraction = %g \n',channels.vol/domainVol);
%{
% create_channels.m were used at the early stages of the code when the
channel input files did not exist
% some of the examples still remain the script and have not been converted
into channel input files
[channels,mesh.heatSourceFunc,soln,u,uL2] = create_channels(mesh.boundary.xi,...
                                                           mesh.boundary.xf,...
                                                            mesh.boundary.yi,...
                                                            mesh.boundary.yf);
%}
%% Analytical distributed heat source
% Note: Only works for NURBS IGFEM !
% To specify analytical distributed heat source for polynomial IGFEM,
% modify the file body_source_functions.cpp in the directory mx_FEM
% mesh.heatSourceFunc = @(x,y) distributed_heat_source(x,y);
mesh.heatSourceFunc = [];
%% analytical or reference solution
isAnalytical = true; % is analytical solution through
%u = @(x,y) analytical_soln(x,y); % for specifying Dirichlet BC based on analytical soln
%soln = @(x,y) analytical_soln(x,y); % for plotting and error analysis

%{
% reference solution based on FEM solution obtained using very fine mesh
isAnalytical = false;
soln = @(x,y) interpolate_soln(meshFine.node.coords,...
      meshFine.elem,UURFine, x,y);

%}            

%% boundary conditions
% 1 : x = xi, 2 : x = xf, 3 : y = yi, 4 : x = yf
mesh.BCs.boundaries = []; % totally insulated boundaries
%mesh.BCs.boundaries = [2,3,4];
%mesh.BCs.types = [1,1,1];
%mesh.BCs.values_or_funcs = {50,50,50};  
isConformingMesh = false;

% specify convective heat coefficient 
%mesh.convect.coef = 8.5+5.71; old value
%mesh.convect.coef = 8+5.6;
%mesh.convect.Tref = 21;
mesh.convect.coef = 0.0; 
mesh.convect.Tref = 0.0;

% for nodes_curves_intersect.m and move_nodes.m 
tol.node = 1e-6; % tolerance for checking if original nodes coincide with channels
tol.boundary = 1e-13; % also for set_boundary_conditions.m
moveNode.distFrac = 0.05; % distance to move node when it coincides with channel
                         % or when a branching point or kink coincides
                         % with an element edge, as a fraction of the
                         % minimum edge length
moveNode.maxAttempts = 5;
moveNode.randDirection = false;

% tolerances for edges_curves_intersect.m
tol.nurbsParam = 1.e-8;
tol.epsco = 1e-15; % computational tolerance for curve_edge_intersect.mexw64
tol.epsge = 1e-6; % geometrical tolerance for curve_edge_intersect.mexw64
tol.cosAngleTol = cos(20*pi/180.0); % cosine of the min angle of the resulting element after edge flipping
tol.halfLineWidthFrac = 1e-4; % half-width fraction of interface
tol.vert = tan(89.9*pi/180); % tan of max positive angle wrt horizontal axis beyond which a line is considered vertical
opt.maxRefineLevel = 0;
opt.refineJuncElem = false; % force refinement of element with branching or kinks
tol.intersectEdges = 1e-13; % for intersect_edges in single_edge_curves_intersect.m


% tolerenace for approximating slender child element with triangle or quadrilateral
% slenderTol = []; % if approximation not desired
slenderTol.minAngle = 5; % min angle of a triangular child element below which NURBS is approximated by line segments
slenderTol.maxAspectRatio = inf;  % max aspect ratio of a quadrialteral child element below which NURBS is approximated by line segments

calcItrsectVel = true; % calculate intersection point velocity

% integration schemes
triNpt1D = 4; % for line integration in regular FEM,poly IGFEM and NURBS IGFEM
triNpt2D = 7; % for element integration in regular FEM and poly IGFEM
quaNpt1Dt = 4; % number of gauss points per knot span along channel for NURBS IGFEM 
quaNpt1Dn = 4; % number of gauss points per knot span normal to channel for NURBS IGFEM
quadNpt1D = 4; % number of gauss points in one direction for quad child element in polynomial IGFEM

% perform FEM analysis
performFEM = true;
polyIGFEM = false;
supg = true; % apply SUPG. Note: do not apply SUPG for constant heat flux model
% plot results
postProcessing = false;
% output paraview file
outfile = 'test';
scalarname = 'T';
% perform error analysismat
errorAnalysis = false;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Start of calculations
%stats = mesh_statistics(mesh.elem.elem_node,mesh.node.coords);
% label edges for refinement later
mesh.elem.elem_node = label(mesh.node.coords,mesh.elem.elem_node); 
%


% the edge_node information is generated by generate_conforming_mesh
fprintf('\ngenerate edge_node\n')
%tic
mesh.edge.edge_node = find_edge_node(mesh.elem.elem_node);
%toc

mesh.edge.length = find_edge_length(mesh.edge.edge_node,mesh.node.coords);

mesh.edge.minLength = min(mesh.edge.length);
tol.halfLineWidth = mesh.edge.minLength*tol.halfLineWidthFrac;
moveNode.dist = mesh.edge.minLength*moveNode.distFrac;

figure
options.showMesh = true;
options.showCtrlPts = false;
options.showDiamVals = false;
options.showNodeNums = false;
plot_mesh_curve(mesh.node.coords, ...
                mesh.elem.elem_node, ...
                channels,options)
% find intersection points and enrichment nodes
if(~isConformingMesh)
    %{
    if (calcItrsectVel)
        fprintf('checking for original nodes that are intersection points \n')
        fprintf('and kinks or branch points on element edges \n')
        %nodeTimer = tic;
        mesh.node.coords ...
                = eliminate_original_node_intersections(mesh.node.coords,...
                                                        mesh.edge.edge_node,...                      
                                                        channels.nurbs,...
                                                        channels.branch_kinks.XX',...
                                                        mesh.boundary,...
                                                        tol.boundary,...
                                                        tol.node,...
                                                        moveNode.maxAttempts,...
                                                        moveNode.dist,...
                                                        moveNode.randDirection);
  
        %toc(nodeTimer);
    end
    %}
    fprintf('\nfinding intersection points \n')
    [mesh.edge,...
     mesh.node,...
     mesh.elem,...
     channels.itrsectParams,...
     refineLevels]...
         =edges_curves_intersect(mesh.edge,mesh.node,...
                                 mesh.elem, ...
                                 channels, ...
                                 tol, ...
                                 refine, ...
                                 otherFlags.calcItrsectVel);
    if refineLevels || isempty(mesh.DT)
        % reconstruct DT if refinement has been carried out
        % or construct DT if DT is empty
        % update other members of mesh.elem
        [mesh.DT,mesh.elem] ...
            = update_mesh_DT_n_elem(mesh.node.coords(1:mesh.node.nOriginalNode,:), ...
                                    mesh.edge.edge_node, ...
                                    mesh.elem);
        
    end
    mesh.elem.elem_edge = find_elem_edge(mesh.elem.elem_node, ...
                                          mesh.edge.edge_node, ...
                                          mesh.node.nOriginalNode);                         
    [mesh.elem.branch_kinks, ...
     mesh.node.coords, ...
     mesh.node.n_node, ...
     mesh.edge.itrsect] ...
                = elem_branching_n_kinks(mesh.DT,...
                                         mesh.elem.elem_edge,...
                                         mesh.edge.itrsect,...
                                         mesh.node.coords,...
                                         mesh.node.n_node,...
                                         channels.branch_kinks,...
                                         channels.itrsectParams,...
                                         tol.nurbsParam, ...
                                         tol.bary, ...
                                         channels.designParamNum, ...
                                         otherFlags.calcItrsectVel);
else
    fprintf('\nconforming mesh\n')

    % find edges sharing a node
    fprintf('generate node_edges\n')
    %tic
    mesh.node.node_edges = find_node_edges(mesh.edge.edge_node,...
                                          size(mesh.node.coords,1));
    %toc
    fprintf('\nfinding intersection points\n')
    %intersectTimer = tic;
    mesh.edge.itrsect = edge_intersect_conforming(mesh.edge.edge_node,...
                                mesh.node.coords, mesh.node.node_edges, ...
                                mesh.elem.elem_edge, mesh.elem.region, channels);
    %toc(intersectTimer);
    
    mesh.elem.junc = [];  
    
    fprintf('\ngenerate elem_edge\n')
    %tic
    mesh.elem.elem_edge = find_elem_edge(mesh.elem.elem_node, ...
                                     mesh.edge.edge_node, ...
                                     size(mesh.node.coords,1));
    %toc
    
end



%
%plot_mesh_curve_itrsect_junc(mesh.node,mesh.elem.elem_node,channels,...
%                   mesh.edge.itrsect,mesh.elem.junc,mesh.elem.kinks,false,false);
%plot_mesh_labels(mesh.node.coords,mesh.elem.elem_node,true,true,mesh.edge.edge_node)
%
%% Set boundary conditions
fprintf('\nsetting boundary conditions \n')
%tic
[mesh.elem, mesh.node] = set_boundary_conditions(mesh.BCs,...
                                                 mesh.elem,...
                                                 mesh.node,...
                                                 mesh.boundary,...
                                                 tol.boundary);
%toc

%% enrichment functions and integration subdomains
% creat parent elements
fprintf('\nconstructing parent elements\n')
%tic
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
                                         mesh.elem.branch_kinks,...
                                         mesh.node.coords,...
                                         mesh.node.nOriginalNode,...
                                         mesh.node.constraint,...
                                         channels,...
                                         tol.nurbsParam,...
                                         polyIGFEM,...
                                         calcItrsectVel);


%toc
%{
elem2plot = [];
plot_mesh_nurbs_parent(mesh.node.coords,mesh.node.nOriginalNode,...
                       mesh.elem.elem_node,mesh.elem.parent,channels,...
                       false,false,elem2plot);
%}

if(mesh.node.n_node==mesh.node.nOriginalNode)
    fprintf('\nconforming mesh\n')
else
    fprintf('\nnon-conforming mesh\n')
    fprintf('constructing child elements\n')
    %tic
    
    [mesh.elem.parent,mesh.node]...edge, Dirichlet, nodeCoords, boundary
        =child_elements_nurbs(mesh.elem.parent,... 
                              mesh.node,...
                              slenderTol,...
                              polyIGFEM);
    %toc                      
    % the global equation number of additional equations for 
    % Lagrange multiplier method                     
    for i = 1:numel(mesh.elem.cstrElems)
        mesh.elem.parent(mesh.elem.cstrElems(i)).cstrRows...
                =  mesh.elem.parent(mesh.elem.cstrElems(i)).cstrRows ...
                 + size(mesh.node.coords,1); 
    end
   
  
    %plot_mesh_nurbs_child(mesh.node.coords,mesh.elem.elem_node,...
    %                      mesh.elem.parent,false,false,elem2plot)
end


%%  Finite Element Code  
%%%%%%%%%%%%%%%%%%%%%%%%%  
% gauss quadrature schemes
%
if(performFEM)
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
    
    %timer = tic;
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
        %                             
        % matlab poly IGFEM    
        %{
        [KPP,KPF,KFP,KFF,PP,PF,UP] = ...
            assemble (mesh, eq_num, n_dof,channels,gauss,supg);
        %}
    else
        gauss.qua1Dt = gauss_points_and_weights(false,quaNpt1Dt,1,'separate'); % along channel
        gauss.qua1Dn = gauss_points_and_weights(false,quaNpt1Dn,1,'separate'); % normal to channel
        fprintf('assembly \n')
        [KPP,KPF,KFP,KFF,PP,PF,UP] = ...
            assemble (mesh, eq_num, n_dof,channels,gauss,supg);
    end
    %toc(timer)
    
    
    
    fprintf('\nsolving equation\n')
    %tic
    [UUR, PUR] = solve_matrix_eqn(KPP,KPF,KFP,KFF,PP,PF,UP, eq_num);  
    %toc            
    
    % output paraview file
    fprintf('\nupdate enrichment node values\n')
    %tic
    UUR2 = update_enrichment_node_value(UUR,mesh.node,mesh.elem,mesh.edge);
    %toc
    
    fprintf('\n total simulation time %g \n',toc(totTimer))
    fprintf('\noutput paraview file\n')
    matlab2vtk_scalar([outfile,'.vtk'],scalarname,mesh.elem,mesh.node,UUR2);
    save(outfile,'mesh','channels','UUR')
    
    %fprintf('\n average temp = %g \n',average_temp(mesh.elem,mesh.node.coords,UUR2));
    % The gauss quadrature for postprocessing (gausspp) can be different 
    % than that for analysis (gauss)
    % Feel free to change it when needed
    gausspp = gauss;
    Tmax = max(UUR2);
    fprintf('\n max temp = %4.12g  \n',Tmax);
    normp = 8;
    gausspp.elem = gauss_points_and_weights(true,16 ,2,'combined');
    pnormT = field_p_norm(mesh.elem,mesh.node.coords,UUR,...
                          normp,Tmax,gausspp);

    fprintf('\n %i-norm temp = %4.12g \n',normp,pnormT);
    [Tvar,Tave,totArea] = field_variance(mesh.elem,mesh.node.coords,UUR,gausspp);
    %Tave = average_temp(mesh.elem,mesh.node.coords,UUR2);
    fprintf('\n Tave = %4.12g , SD = %4.12g  \n',Tave,sqrt(Tvar));
   
    nu = kinematic_viscosity(channels.viscosity,channels.density,Tave, ...
                             channels.viscosityModel);
    [pressure,mass] = network_pressure_mass_flow_rate(channels.contvty,...
                                                      channels.nurbs,...
                                                      channels.diams,...
                                                      channels.heights,...
                                                      nu,...
                                                      channels.inletEndPoint,...
                                                      channels.massin,...
                                                      channels.powerXdensity,...
                                                      channels.pressureOutletEndPoint,...
                                                      channels.pressureOut,...
                                                      channels.crossSection);
    fprintf('\n Inlet pressure(s) = %g\n',pressure(channels.inletEndPoint));                                              
    %figure                                              
    %plot_channel_network(channels.contvty,channels.nurbs,pressure/1000.0,mass*1e6*60/1065)       
    %% postprocessing
    if(postProcessing)    
        scale.length = 1.0;
        scale.uo = 0.0;
        scale.ud = 1.0;
        fprintf('\npost processing\n')
        % plot solution    
        specs.line = 'b-';
        specs.marker = 'bo';
        specs.mksize = 6.0;
        specs.analytical = 'r-.';
        specs.width = 2;
        npts =1; % number of points between pairs of intersections between the line and the mesh
        if (~isempty(soln))
            line = [0.0,0.01,1,0]; % line = [x0,y0,dx,dy]   
            figure
            plot_interpolated_n_analytical_soln_along_line(mesh.node.coords,...
                mesh.elem,mesh.edge,UUR,line,npts,tol.epsge,specs,isAnalytical,soln);
        else  
            
            line = [0.0,0.18,1,0];
            figure
            plot_interpolated_n_analytical_soln_along_line(mesh.node.coords,...
                mesh.elem,mesh.edge,UUR,line,npts,tol.epsge,specs);      


        end
        %plot_fill_elements_node_wise_quantity(mesh.node.coords,mesh.elem,UUR2,scale)
    end
end
%
%% error analysis

if(errorAnalysis)
    fprintf('\nperforming error analysis\n')
    [error.totEL2,error.totEH1,error.totParentEL2,error.totParentEH1,...
     error.totL2norm,error.totH1norm,mesh.elem]=error_L2_H1(soln,UUR,mesh.node.coords,mesh.elem,gauss);
    error.numberParent =  nnz(cat(1,mesh.elem.parent(:).type) > 1);
    error.maxEdgeLength = max(mesh.edge.length);
    error.minEdgeLength = min(mesh.edge.length);
    error.aveEdgeLength = mean(mesh.edge.length);
    error.n_node = mesh.node.n_node;
    n_pre_temp = length(unique(mesh.node.Dirichlet.temp_node));
    if(n_pre_temp ~= mesh.node.Dirichlet.n_pre_temp)
        warning('some Dirichlet nodes are not unique')
    end
    error.dof = mesh.node.n_node - n_pre_temp;    
    error.normalizedTotEL2 = error.totEL2/error.totL2norm;            
    uH1 = 1.0;
    %error.normalizedTotEH1 = error.totEH1/error.totH1norm;
    error.normalizedTotEH1 = error.totEH1/uH1;
    error
    %plot_fill_elements_elem_wise_error_or_norm(mesh.node.coords,mesh.elem,'L2','error');
end

fclose ('all');

toc(totTimer)