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
path(path, '../nurbs_toolbox')
path(path, './mesh_conforming_abaqus')
path(path, './mesh')
path(path, '../SISL/mx_SISL') % for windows
path(path, './ChannelFiles')
path(path, './InputFiles')
path(path, '../microchannels-optimizer-2D/M_optimization')
path(path, './export_fig')
%path(path, '../../SISL_linux') % for LINUX

totTimer = tic;
%% MESH AND USER INPUT
inputFile = 'cubesat_panel.in';
 [mesh,gauss,tol,refine,...
  otherFlags,postprocess, ...
  moveNode] = read_inputs(inputFile);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Other user inputs that must be specified regardless of choice I or II
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Interface, analytical body source (if available) and analytical 
% solution (if available)
%channelFile = '3x3localRefLp075deg6.channel';
%channelFile = '3x3Lp075deg6_nBlk1_optimal.channel';
channelFile = 'cubesat_bifurcating4.channel';
options.suppressDesignParams = false;
options.selfIntersectTol = tol.channelSelfIntersect;
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

% perform error analysismat
errorAnalysis = false;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Start of calculations
%stats = mesh_statistics(mesh.elem.elem_node,mesh.node.coords);
% label edges for refinement later
% IMPORTANT: This causes elem_node to be inconsistent with
%            DT.ConnectivityList
if refine.maxRefineLevel
    mesh.elem.elem_node = label(mesh.node.coords,mesh.elem.elem_node); 
end

tol.halfLineWidth = mesh.edge.minLength*tol.halfLineWidthFrac;
moveNode.dist = mesh.edge.minLength*moveNode.distFrac;


            
%plot_mesh_labels(mesh.node.coords,mesh.elem.elem_node, ...
%                 false,false,mesh.edge.edge_node)
% find intersection points and enrichment nodes
if(~otherFlags.isConformingMesh)
    nodesMoved = false;
    %{
    if otherFlags.calcItrsectVel
        fprintf('checking for original nodes that are intersection points \n')
        fprintf('and kinks or branch points on element edges \n')
        %nodeTimer = tic;
        [mesh.node.coords,nodesMoved] ...
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
  
        %IMPORTANT: at this point, if mesh.node.coords is changed, it will
        % no longer be consistent with mesh.DT.Points. 
        % So mesh.DT should be updated before being used
        % However, if refinement is requested, mesh.elem.elem_node will be
        % ordered in a certain way for the bisection method.
        % So before the refinement, do not change mesh.elem.elem_node
    end
    %}
    fprintf('\nfinding intersection points \n')
    %tic    
    [mesh.edge,...
     mesh.node,...
     mesh.elem,...
     channels.itrsectParams,...
     refineLevel]...
         =edges_curves_intersect(mesh.edge,mesh.node,...
                                 mesh.elem, ...
                                 channels, ...
                                 tol, ...
                                 refine, ...
                                 otherFlags.calcItrsectVel);  
                             
    if refineLevel || nodesMoved
        fprintf('Update DT\n')
        % reconstruct DT if refinement has been carried out or
        % the nodes have been moved
        [mesh.DT,mesh.elem] ...
            = update_mesh_DT_n_elem(mesh.node.coords(1:mesh.node.nOriginalNode,:), ...
                                    mesh.edge.edge_node, ...
                                    mesh.elem);
        mesh.elem.elem_edge = find_elem_edge(mesh.elem.elem_node, ...
                                             mesh.edge.edge_node, ...
                                             mesh.node.nOriginalNode);  
    elseif refine.maxRefineLevel
        % although no refinement has been carried out, elem_node has 
        % become inconsistent with DT after relabeling by the label
        % function for bisection. So, we need to 
        % revert it back to its original connectivity
        mesh.elem.elem_node = mesh.DT.ConnectivityList;
    end
    
    figure
    options.showMesh = true;
    options.showCtrlPts = false;
    options.showDiamVals = false;
    options.showNodeNums = false;
    options.showElemNums = false;
    options.showElemVals = false;
    plot_mesh_curve(mesh.node.coords, ...
                    mesh.elem.elem_node, ...
                    mesh.elem.heatSource, ...
                    channels,options)
            
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
    % nOriginalEnrichNode only includes enrichment nodes that arise due to 
    % the intersections
    % If NURBS-IGFEM is used, additional enrichment nodes corresponding to
    % additional control points may arise
    mesh.node.nOriginalEnrichNode = mesh.node.n_node;                                 
    %toc
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
                                mesh.node.coords, mesh.node.node_edge, ...
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
                                         otherFlags.calcItrsectVel);
%
[mesh.elem.parent, ...
 mesh.elem.cstrElems, ...
 mesh.elem.nIGFEMelems] ...
                 = parent_elements_nurbs(mesh.elem.parent,...
                                         mesh.elem.elem_node,...
                                         mesh.elem.material,...
                                         mesh.material,...
                                         mesh.elem.branch_kinks,...
                                         mesh.node.coords,...
                                         mesh.node.nOriginalNode,...
                                         mesh.node.constraint,...
                                         channels,...
                                         tol.nurbsParam,...
                                         otherFlags.polyIGFEM,...
                                         otherFlags.calcItrsectVel);


%toc
%
elem2plot = [];
plot_mesh_nurbs_parent(mesh.node.coords,mesh.node.nOriginalNode,...
                       mesh.elem.elem_node,mesh.elem.parent,channels,...
                       false,false,elem2plot);
%    
if(mesh.node.n_node==mesh.node.nOriginalNode)
    fprintf('\nconforming mesh\n')
else
    fprintf('\nnon-conforming mesh\n')
    fprintf('constructing child elements\n')
    %tic
    
    [mesh.elem.parent,mesh.node]...edge, Dirichlet, nodeCoords, boundary
        =child_elements_nurbs(mesh.elem.parent,... 
                              mesh.node,...
                              tol.collinear,...
                              tol.slender,...                            
                              otherFlags.polyIGFEM);
    %toc                      
    % the global equation number of additional equations for 
    % Lagrange multiplier method                     
    for i = 1:numel(mesh.elem.cstrElems)
        mesh.elem.parent(mesh.elem.cstrElems(i)).cstrRows...
                =  mesh.elem.parent(mesh.elem.cstrElems(i)).cstrRows ...
                 + size(mesh.node.coords,1); 
    end
    
    %{
    elem2plot = [];
    plot_mesh_nurbs_child(mesh.node.coords,mesh.elem.elem_node,...
                          mesh.elem.parent,true,true,elem2plot)
    %}
end


%%  Finite Element Code  
%%%%%%%%%%%%%%%%%%%%%%%%%  
% gauss quadrature schemes
%
if(otherFlags.performFEM)
    fprintf('\nperforming FEM\n')

    % eq_num:  Equation number assigned to each node
    % n_dof: The number of degree of freedom in the model
    [n_dof, eq_num] = initialize(size(mesh.node.coords,1) ...
                                 +numel(mesh.node.constraint.temp_node), ...
                                 mesh.node.Dirichlet.n_pre_temp, ...
                                 mesh.node.Dirichlet.temp_node);

    %timer = tic;
    if (otherFlags.polyIGFEM)
        fprintf('\nassembly of polynomial IGFEM stiffness matrix\n')
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
                                         otherFlags.supg);
    else
        fprintf('assembly \n')
        [KPP,KPF,KFP,KFF,PP,PF,UP] = ...
            assemble (mesh, eq_num, n_dof,channels,gauss,otherFlags.supg);
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
    matlab2vtk_scalar([postprocess.outfile,'.vtk'],postprocess.scalarname, ...
                      mesh.elem,mesh.node,UUR2);
    save(postprocess.outfile,'mesh','channels','UUR')
    
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
    if(postprocess.matlab)    
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
