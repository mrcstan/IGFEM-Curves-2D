%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Marcus Tan on 11/16/2013
%%% Last modified date: 11/16/2013
%%% Copyright 2013 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this function defines the vertices of the boundary of the domain
function boundary = define_boundary_vertices(xi,xf,yi,yf)
    boundary.xv = [xi,xf,xf,xi];
    boundary.yv = [yi,yi,yf,yf];
end