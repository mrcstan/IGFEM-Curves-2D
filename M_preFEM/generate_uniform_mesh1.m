%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Marcus Tan on 7/2/2014
%%% Last modified date: 7/2/2014
%%% Copyright 2013 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [elem_nodes, node_coords] = generate_uniform_mesh(m,n,xi,xf,yi,yf,type)
m = m+1;
n = n+1;
[X,Y] = meshgrid(linspace(xi,xf,m),linspace(yi,yf,n));

switch type
    case 1
        X = X';
        Y = Y'; 
        A = repmat(1:m:(m*(n-1)),m-1,1);
        SW = A(:)+kron(ones(n-1,1),(1:m-1)')-1;
        SE = SW + 1;
        NW = SW + m;
        NE = NW + 1;
        elem_nodes = zeros(2*(m-1)*(n-1),3);
        elem_nodes(1:2:end) = [SW,NE,NW];
        elem_nodes(2:2:end) = [SW,SE,NE];
        
    case 2
        elem_nodes = delaunay(X,Y);
    otherwise
        error('unrecognized option')
end
node_coords = [X(:),Y(:)];
end
