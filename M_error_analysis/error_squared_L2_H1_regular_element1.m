%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Marcus Tan on 7/26/2013
%%% Last modified date: 8/13/2013
%%% Copyright 2013 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUT:
%   element-wise quantity
%   varargout(1): (L2 error)^2
%            (2): (L2 norm of analytical solution)^2 
%            (3): (H2 error)^2
%            (4): (H2 norm of analytical solution)^2
function varargout= error_squared_L2_H1_regular_element(analytical,UUR,nodeCoords,...
                                              elem_node,gauss)
% initalize quantities of interest
varargout{1}=0;
varargout{2}=0;
if(nargout>2)
    varargout{3}=0;
    varargout{4}=0;
end

Xel = [nodeCoords(elem_node,1), ...
       nodeCoords(elem_node,2)]';
shape = 1;

for j = 1:gauss.tri2D.npt   
    r=gauss.tri2D.pt(j,1:2);
    w=gauss.tri2D.weight(j);
    [N, DN] = shape_funct(r, shape);
    J = Xel*DN; % Xel is a 2x3 matrix, DN is a 3x2 matrix
    XX=Xel*N;
    UU=UUR(elem_node)'*N;  
    u_du=analytical(XX(1),XX(2));
    if(isempty(u_du))
        warning('integration point cannot be found in any element')
    end
    factor=det(J)*w;
    UU_ua_sq=(UU-u_du(1))^2;
    varargout{1}=varargout{1}+UU_ua_sq*factor;
    ua_sq=u_du(1)^2;
    varargout{2}=varargout{2}+ua_sq*factor;    
    if(nargout>2)
        B=DN/J;
        dU=UUR(elem_node)'*B;   
        vector=[dU(1)-u_du(2),dU(2)-u_du(3)];
        varargout{3}=varargout{3}+(UU_ua_sq+dot(vector,vector))*factor;
        vector=[u_du(2),u_du(3)];
        varargout{4}=varargout{4}+(ua_sq+dot(vector,vector))*factor;
    end
end
end

