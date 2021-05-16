%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Marcus Tan on 9/1/2013
%%% Last modified date: 09/1/2013
%%% Copyright 2013 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function outputs the 1D shape functions and its derivatives
function [N, DN] = shape_funct_1D(r)
    N=[0.5*(1-r);0.5*(1+r)];
    DN=[-0.5;0.5];
end