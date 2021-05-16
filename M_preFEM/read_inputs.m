%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Marcus Tan on 3/17/2016
%%% Copyright 2016 University of Illinois at Urbana-Champaign. 
%%% All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read channel file input
% INPUT: 
%   inputFile: input file name 
%   options: the options specified here, if non-empty, take precedence 
%            over those in the input file
%       options.meshSizes: number of elements in along the axes
%       options.outfile: vtk and MATLAB variable file name
% OUTPUT:
%   the outputs are quite self-explanatory
% FORMAT:
%{
A line starting with # will be ignored.
Each keyword is to be followed by a comma and then value(s).
The general structure is better explained with an example file as shown below.
Most keywords are explained in the comments of the example file.
Here are a few keywords that need more explanation:
mesh_orientation, <value>
   0: Delaunay triangulation
   1: diagonal along NW direction
   2: diagonal along NE direction
BC_boundaries, <values>
   <values> consists of a combination of the following:
   1 : x = xi, 2 : x = xf, 3 : y = yi, 4 : y = yf
   insulated: % totally insulated boundaries
BC_types, <values>
   <values> consists of a combination of the following:
   1: Dirichlet
   2: Neumann/convection
   if BC_boundaries, insulated, just put any value
BC_heat_flux_or_amb_temp, <values>
   <values> consists of a combination of the following: 
   If type == 2 and Neumann required, values_or_funcs = q''
   if BC_boundaries, insulated, just put any value
%}
% REMARKS: i) This function has not been tested with an Abaqus mesh file
%          ii) mesh.BCs.values_or_funcs is a cell as it can also accept 
%               function handles (this is not available in the input file)
%--------------------------------------------------------------------------
% Example input files can be found in the directory InputFiles.
% Here is one example.
% EXAMPLE 1:
%{
# if meshFile = NA, use built in mesh generator
# otherwise, assume it is an Abaqus mesh file
mesh_file, NA
# domain_bounds, xi xf yi yf
domain_bounds, 0 0.15 0 0.2
# mesh_sizes, nx ny nz (options.meshSizes would take precedence over these)
mesh_sizes, 15 20
# mesh_orientation (see read_inputs.m)
mesh_orientation, 0
# material conductivity
conductivity, 2.7
thickness, 0.003
# heat source 
heat_source, 500
# convection_coefficient, 13.6
convection_coefficient, 0.0
ambient_temperature, 21.0
# boundary conditions (see read_input.m)
BC_boundaries, insulated
BC_types, insulated
BC_values, insulated
# tolerance for NURBS parameters
tol_nurbs_parameter, 1e-8
# computational tolerance for curve_edge_intersect
tol_computational, 1e-15
# geometrical tolerance for curve_edge_intersect
tol_geometrical, 1e-5
# half-width fraction of channel
tol_half_line_width_fraction, 1e-4
# tolerance for intersect_edges in single_edge_curves_intersect.m
tol_intersect_edges, 1e-10
# min angle of the resulting elements after edge flipping
tol_edge_flip_angle, 20
# max positive angle wrt horizontal axis beyond which a line
# is considered vertical
tol_vertical_angle, 89.9
# min angle of a triangular child element below which NURBS is approximated by line segments
tol_slender_min_angle, 5
# max aspect ratio of a quadrialteral child element below which NURBS is approximated by line segments
tol_slender_max_aspect_ratio, inf 
# tolerance for checking if original nodes coincide with channels
tol_node, 1e-6
# tolerance for determining whether nodes lie on domain boundaries
tol_boundary, 1e-13
# tolerance for determining if channels self-intersect
tol_self_intersect, 1e-10
# barycentric coordinate tolerance for finding branching points and kinks
tol_barycentric, 1e-10
# maximum number of refinements
refine_max_levels, 0
# force refinement of elements containing branching points
refine_branch_elem, false
# distance to move node when it coincides with channel
# or when a branching point or kink coincides
# with an element edge, as a fraction of the
# minimum edge length
move_node_dist_frac, 0.05
# number of attempts to move node
move_node_max_attempts, 10
# randomize directions in which the nodes are moved
move_node_rand_direction, true
# gauss quadrature for line integration in regular, poly and NURBS IGFEM
gauss_line_num_pts, 4
# gauss quadrature for regular elements
gauss_regular_elem_num_pts, 7
# number of gauss points per knot span along channel for NURBS IGFEM 
gauss_along_channel, 4
# number of gauss points per knot span normal to channel for NURBS IGFEM
gauss_normal_to_channel, 4
# number of gauss points in one direction for quad child element in polynomial IGFEM (disabled)
gauss_quad_each_direction, 4
# other flags
is_conforming_mesh, false
calculate_intersection_velocities, true
apply_supg, true
perform_FEM, true
# performs polynomial IGFEM instead of NURBS IGFEM
polynomial_IGFEM, true
# output vtk file name without extension (options.outfile takes precedences over this)
vtk_file_name, test
# output vtk scalar name
vtk_scalar_name, T_C
# postprocess in MATLAB
postprocess, false
%}
function [mesh,gauss,tol,refine, ...
          otherFlags,postprocess,...
          moveNode] = read_inputs(inputFile,options)                                   

if (~ischar(inputFile))
    error('input file name must be a character array');
end

if nargin < 2
    options = struct;
end

if ~isfield(options,'meshSizes')
    options.meshSizes = [];
else
    if numel(options.meshSizes) ~= 3
        error('options.meshSizes must be a vector of three entries')
    end
end

if ~isfield(options,'outfile')
    options.outfile = [];
end
fileID = fopen(inputFile,'r');
line = fgetl(fileID);

% first keyword must be mesh_file
meshFile = nan;
while ischar(line)
    line = strtrim(line); 
    % skip comment line (starting with #) and empty line
    if (isempty(line) ||strncmpi(line,'#',1))
        line = fgetl(fileID);
        continue
    end
    splitStr = regexp(line,',','split');
    splitStr{1} = strtrim(splitStr{1});
    if numel(splitStr{2}) < 2
        error('keyword should be followed by a comma and values or string')
    end
    splitStr{2} = strtrim(splitStr{2});
    if strcmpi(splitStr{1},'mesh_file')       
        if strcmp(splitStr{2},'NA')
            meshFile = [];
        elseif ischar(splitStr{2})
            meshFile = splitStr{2};
        else
            error('a string should followed mesh_file,')
        end
        break
    end
    line = fgetl(fileID);
end

if isnan(meshFile)
    error('mesh_file keyword not found')
elseif ischar(meshFile)
    mesh = read_abaqus_mesh(meshFile);
end

while ischar(line)
    line = strtrim(line);
    
    % skip comment line (starting with #) and empty line
    if (isempty(line) ||strncmpi(line,'#',1))
        line = fgetl(fileID);
        continue
    end
    splitStr = regexp(line,',','split');
    splitStr{1} = strtrim(splitStr{1});
    if numel(splitStr{2}) < 2
        error('keyword should be followed by a comma and values or string')
    end    
    splitStr{2} = strtrim(splitStr{2});
    if strcmpi(splitStr{1},'domain_bounds')
        domainBounds = str2num(splitStr{2});
        if numel(domainBounds) ~= 4
            error('domain_bounds keyword should be followed by 4 values')
        end
        if (domainBounds(2) <= domainBounds(1) ...
            || domainBounds(4) <= domainBounds(3))
            error('some entries of domainBounds are not in expected order')
        end
    elseif strcmpi(splitStr{1},'mesh_sizes')
        meshSizes = str2num(splitStr{2});
        if numel(meshSizes) ~= 2
            error('mesh_sizes keyword should be followed by 2 values')
        end
        if any(meshSizes < 1)
            error('mesh_sizes all entries must be > 1')
        end
    elseif strcmpi(splitStr{1},'mesh_orientation')
        meshOrientation = str2double(splitStr{2});
    elseif strcmpi(splitStr{1},'conductivity')
        conductivity = str2double(splitStr{2});
    elseif strcmpi(splitStr{1},'thickness')
        thickness = str2double(splitStr{2});    
    elseif strcmpi(splitStr{1},'heat_source')
        heatSource_n_bounds = str2num(splitStr{2});
        if numel(heatSource_n_bounds) ~= 1 ...
          && numel(heatSource_n_bounds) ~= 5      
            error('heat_source keyword must be followed by either 1 or 5 values')
        end
    elseif strcmpi(splitStr{1},'convection_coefficient')
        convecCoef = str2double(splitStr{2});
    elseif strcmpi(splitStr{1},'ambient_temperature')
        convecTamb = str2double(splitStr{2});    
    elseif strcmpi(splitStr{1},'BC_boundaries')
        if strcmpi(splitStr{2},'insulated')
            BC_boundaries = [];
        else
            BC_boundaries = str2num(splitStr{2});
        end
    elseif strcmpi(splitStr{1},'BC_types')
        BC_types = str2num(splitStr{2});
    elseif strcmpi(splitStr{1},'BC_values')
        BC_values = str2num(splitStr{2});
    elseif strcmpi(splitStr{1},'tol_nurbs_parameter')
        tol.nurbsParam = str2double(splitStr{2});  
    elseif strcmpi(splitStr{1},'tol_computational')
        tol.epsco = str2double(splitStr{2});  
    elseif strcmpi(splitStr{1},'tol_geometrical')
        tol.epsge = str2double(splitStr{2});    
    elseif strcmpi(splitStr{1},'tol_half_line_width_fraction')
        tol.halfLineWidthFrac = str2double(splitStr{2});
    elseif strcmpi(splitStr{1},'tol_intersect_edges')
        tol.intersectEdges = str2double(splitStr{2});    
    elseif strcmpi(splitStr{1},'tol_edge_flip_angle')
        tol.cosAngleTol = cos(str2double(splitStr{2})*pi/180);
    elseif strcmpi(splitStr{1},'tol_vertical_angle')
        tol.vert = str2double(splitStr{2});
    elseif strcmpi(splitStr{1},'tol_collinear')
        tol.collinear = str2double(splitStr{2});    
    elseif strcmpi(splitStr{1},'tol_slender_min_angle')
        tol.slender.minAngle = str2double(splitStr{2});
    elseif strcmpi(splitStr{1},'tol_slender_max_aspect_ratio')
        tol.slender.maxAspectRatio = str2double(splitStr{2});    
    elseif strcmpi(splitStr{1},'tol_node')
        tol.node = str2double(splitStr{2});
    elseif strcmpi(splitStr{1},'tol_boundary')
        tol.boundary = str2double(splitStr{2});
    elseif strcmpi(splitStr{1},'tol_self_intersect')
        tol.channelSelfIntersect = str2double(splitStr{2});
    elseif strcmpi(splitStr{1},'tol_barycentric')
        tol.bary = str2double(splitStr{2});    
    elseif strcmpi(splitStr{1},'refine_max_levels')
        refine.maxRefineLevel = str2double(splitStr{2});  
    elseif strcmpi(splitStr{1},'refine_branch_elem')
        refine.refineJuncElem = str2logical(splitStr{2});     
    elseif strcmpi(splitStr{1},'move_node_dist_frac')
        moveNode.distFrac = str2double(splitStr{2});
    elseif strcmpi(splitStr{1},'move_node_max_attempts')
        moveNode.maxAttempts = str2double(splitStr{2});
    elseif strcmpi(splitStr{1},'move_node_rand_direction')
        moveNode.randDirection = str2logical(splitStr{2});  
    elseif strcmpi(splitStr{1},'gauss_line_num_pts')
        triNpt1D = str2double(splitStr{2});
    elseif strcmpi(splitStr{1},'gauss_regular_elem_num_pts')
        triNpt2D = str2double(splitStr{2});
    elseif strcmpi(splitStr{1},'gauss_along_channel')
        quaNpt1Dt = str2double(splitStr{2});
    elseif strcmpi(splitStr{1},'gauss_normal_to_channel')
        quaNpt1Dn = str2double(splitStr{2});  
    elseif strcmpi(splitStr{1},'gauss_quad_each_direction')
        quadNpt1D = str2double(splitStr{2});        
    elseif strcmpi(splitStr{1},'is_conforming_mesh')
        otherFlags.isConformingMesh = str2logical(splitStr{2});
    elseif strcmpi(splitStr{1},'calculate_intersection_velocities')
        otherFlags.calcItrsectVel = str2logical(splitStr{2});
    elseif strcmpi(splitStr{1},'apply_supg')
        otherFlags.supg = str2logical(splitStr{2});
    elseif strcmpi(splitStr{1},'perform_FEM')
        otherFlags.performFEM = str2logical(splitStr{2});
    elseif strcmpi(splitStr{1},'polynomial_IGFEM')
        otherFlags.polyIGFEM = str2logical(splitStr{2});    
    elseif strcmpi(splitStr{1},'vtk_file_name')
        if ~ischar(splitStr{2})
            error('value after vtk_file_name, should be a string')
        end
        postprocess.outfile = splitStr{2};    
    elseif strcmpi(splitStr{1},'vtk_scalar_name')
        if ~ischar(splitStr{2})
            error('value after vtk_scalar_name, should be a string')
        end
        postprocess.scalarname = splitStr{2};
    elseif strcmpi(splitStr{1},'postprocess')  
        postprocess.matlab = str2logical(splitStr{2});
    end
    line = fgetl(fileID);
end

if isempty(meshFile)
    % the mesh sizes specified as part of the options in 
    % a MATLAB script takes precedence over those specified in the input
    % file
    if ~isempty(options.meshSizes)
        meshSizes = options.meshSizes;
    end
    % built in structured mesh generator
    [mesh.elem.elem_node, ...
     mesh.node.coords, ...
     mesh.edge.edge_node, ...
     mesh.DT] ...
            = generate_uniform_mesh(meshSizes,... 
                                    domainBounds,...
                                    meshOrientation); 
    mesh.boundary.xi = domainBounds(1);
    mesh.boundary.xf = domainBounds(2);
    mesh.boundary.yi = domainBounds(3);
    mesh.boundary.yf = domainBounds(4);
    mesh.meshSizes = meshSizes;
elseif ischar(meshFile)
    % Abaqus mesh file
    warning('read from Abaqus mesh file not yet checked')
    mesh = read_abaqus_mesh(meshFile);
    mesh.edge.edge_node = find_edge_node(mesh.elem.elem_node);
    mesh.DT = delaunayTriangulation(mesh.node.coords, ...
                                    mesh.edge.edge_node);
    mesh.elem.elem_node = mesh.DT.ConnectivityList;
    mesh.boundary.xi = min(mesh.node.coords(:,1));
    mesh.boundary.xf = max(mesh.node.coords(:,1));
    mesh.boundary.yi = min(mesh.node.coords(:,2));
    mesh.boundary.yf = max(mesh.node.coords(:,2));
end
mesh.edge.length = find_edge_length(mesh.edge.edge_node,mesh.node.coords);
mesh.edge.minLength = min(mesh.edge.length);

mesh.elem.elem_edge = find_elem_edge(mesh.elem.elem_node, ...
                                     mesh.edge.edge_node, ...
                                     size(mesh.node.coords,1));
                                 
mesh.domainArea = (mesh.boundary.xf - mesh.boundary.xi) ...
                  *(mesh.boundary.yf - mesh.boundary.yi);
mesh.domainVol =  mesh.domainArea*thickness;   

mesh.material.conductivity = conductivity*thickness; 
mesh.elem.material = int32(ones(size(mesh.elem.elem_node,1),1));
% distributed heat source 
nElems = size(mesh.elem.elem_node,1);
if numel(heatSource_n_bounds) == 1
    mesh.elem.heatSource = heatSource_n_bounds*ones(nElems,1);
elseif numel(heatSource_n_bounds) == 5
    if (heatSource_n_bounds(2) >= heatSource_n_bounds(3) ...
        || heatSource_n_bounds(4) >= heatSource_n_bounds(5))
        error('last 4 values of heat_source should be xi xf yi yf where xi<xf,yi<yf')
    end
    mesh.elem.heatSource = zeros(nElems,1);
    for i = 1:nElems
        Xel = mesh.node.coords(mesh.elem.elem_node(i,:),:);
        if all(Xel(:,1) >= heatSource_n_bounds(2) ...
               & Xel(:,1) <= heatSource_n_bounds(3) ...
               & Xel(:,2) >= heatSource_n_bounds(4) ...
               & Xel(:,2) <= heatSource_n_bounds(5))
           mesh.elem.heatSource(i) = heatSource_n_bounds(1);
        end
    end
end
% convection coefficient and ambient temperature    
mesh.convect.coef = convecCoef; 
mesh.convect.Tref = convecTamb;

nBCs = numel(BC_boundaries);
if numel(BC_types) ~= nBCs
    error('number of BC_types entries must be equal that of BC_boundaries')
elseif numel(BC_values) ~= nBCs
    error('number of BC_heat_flux_or_temperature entries must be equal that of BC_boundaries')
end
mesh.BCs.boundaries = BC_boundaries;
if isempty(BC_boundaries)
    mesh.BCs.types = [];
    mesh.BCs.values_or_funcs = [];
else
    mesh.BCs.types = BC_types;
    mesh.BCs.values_or_funcs = cell(1,nBCs);
    for i = 1:nBCs
        mesh.BCs.values_or_funcs{i} = BC_values(i);

    end
end
gauss.line = gauss_points_and_weights(true,triNpt1D,1,'combined');
gauss.elem = gauss_points_and_weights(true,triNpt2D ,2,'combined');

if otherFlags.polyIGFEM
    gauss.quadElem = gauss_points_and_weights(false,quadNpt1D,2,'combined');
else
    gauss.qua1Dt = gauss_points_and_weights(false,quaNpt1Dt,1,'separate'); % along channel
    gauss.qua1Dn = gauss_points_and_weights(false,quaNpt1Dn,1,'separate'); % normal to channel
end

% options.outfile takes precedence over postprocess.outfile if it is a
% string
if ischar(options.outfile)
    postprocess.outfile = options.outfile;
end

end

function logical = str2logical(str)
if ~ischar(str)
    error('input must be string')
end
if strcmpi(str,'true')
    logical = true;
elseif strcmpi(str,'false')
    logical = false;
else
    error('string must be either true or false')
end
end
