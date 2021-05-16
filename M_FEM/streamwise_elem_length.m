%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Marcus Tan on 7/24/2014
%%% Last modified date: 7/24/2014
%%% Copyright 2014 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this function calculates h_e, the element length along the streamline
% direction for an element
% INPUT: 
%   unitTan: [tx;ty] unit tangent in streamwise direction
%   B: [dN1/dx,dN1/dy;dN2/dx,dN2/dy;...;dNn/dx,dNn/dy];
% OUTPUT:
%   he
%   Bsw: a length n column vector of the derivatives of shape functions in streamwise direction
function [he,Bsw] = streamwise_elem_length(unitTan,B)
Bsw = B*unitTan;
he = 2.0/sum(abs(Bsw));
end