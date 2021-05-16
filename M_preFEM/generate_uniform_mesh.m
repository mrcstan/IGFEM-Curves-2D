%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Marcus Tan on 7/2/2014
%%% Last modified date: 10/31/2015
%%% Copyright 2013 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function generates a structured mesh of a rectangular domain in 2D
% or a cuboidal domain in 3D
% INPUT:
%   nel = [nx,ny] in 2D ([nx,ny,nz] in 3D) of the number of elements in
%         x,y (x,y,z) directions
%   xx = [xi,xf,yi,yf] in 2D ([xi,xf,yi,yf,zi,zf] in 3D) of the bounding
%   coordinates of the rectangular (cuboidal) domain
%   type: 
%           0: Delaunay triangulation
%           1: diagonal along NW direction (NA for 3D)
%           2: diagonal along NE direction (NA for 3D)
function [elem_nodes,node_coords,edge_nodes,DT] = generate_uniform_mesh(nel,xx,type)
nel = nel + 1;
if 2*numel(nel) ~= numel(xx)
    error('input for number of elements must be consistent with the bounding coordinates')
end
% 2D case
if numel(nel) == 2
    [X,Y] = meshgrid(linspace(xx(1),xx(2),nel(1)),linspace(xx(3),xx(4),nel(2)));
    switch type
        case 0
            DT = delaunayTriangulation(X(:),Y(:));
            elem_nodes = DT.ConnectivityList;
            edge_nodes = edges(DT);
            node_coords = [X(:),Y(:)];
        case 1
            %elem_nodes = delaunay(X,Y); % does not guarantee a NW
            %orientation due to floating point round off error
            node_coords = [reshape(X',[],1),reshape(Y',[],1)];
            A = repmat(1:nel(1):(nel(1)*(nel(2)-1)),nel(1)-1,1);
            SW = A(:)+kron(ones(nel(2)-1,1),(1:nel(1)-1)')-1;
            SE = SW + 1;
            NW = SW + nel(1);
            NE = NW + 1;
            elem_nodes = zeros(2*(nel(1)-1)*(nel(2)-1),3);
            elem_nodes(1:2:end) = [SW,SE,NW];
            elem_nodes(2:2:end) = [SE,NE,NW];
            edge_nodes = find_edge_node(elem_nodes);
            DT = delaunayTriangulation(node_coords, ...
                                       edge_nodes);
            elem_nodes = DT.ConnectivityList;
            
        case 2
            node_coords = [reshape(X',[],1),reshape(Y',[],1)];
            A = repmat(1:nel(1):(nel(1)*(nel(2)-1)),nel(1)-1,1);
            SW = A(:)+kron(ones(nel(2)-1,1),(1:nel(1)-1)')-1;
            SE = SW + 1;
            NW = SW + nel(1);
            NE = NW + 1;
            elem_nodes = zeros(2*(nel(1)-1)*(nel(2)-1),3);
            elem_nodes(1:2:end) = [SW,NE,NW];
            elem_nodes(2:2:end) = [SW,SE,NE];
            edge_nodes = find_edge_node(elem_nodes);
            DT = delaunayTriangulation(node_coords, ...
                                       edge_nodes);
            elem_nodes = DT.ConnectivityList;             
        otherwise
            error('unrecognized type of mesh orientation')
    end

   
% 3D case
elseif numel(nel) == 3
    [X,Y,Z] = meshgrid(linspace(xx(1),xx(2),nel(1)), ...
                       linspace(xx(3),xx(4),nel(2)), ...
                       linspace(xx(5),xx(6),nel(3)));
    node_coords = [X(:),Y(:),Z(:)];
    switch type
        case 0 % delaunay triangulation
            %elem_nodes = delaunay(X(:),Y(:),Z(:));
            DT = delaunayTriangulation(X(:),Y(:),Z(:));
            elem_nodes = DT.ConnectivityList;
            edge_nodes = edges(DT);
        otherwise
            error('unrecognized type of mesh orientation')
    end
end
end
