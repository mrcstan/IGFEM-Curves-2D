%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Marcus Tan on 12/26/2014
%%% Copyright 2014 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this function computes the error in a polynomial IGFEM element 
function varargout =compute_polyIGFEM_element_error(refSoln,...
                                                    UUR,...
                                                    nodeCoords,...
                                                    parent,...
                                                    gauss)    
varargout{1} = 0;
varargout{2} = 0;
if (nargout > 2)
    varargout{3} = 0;
    varargout{4} = 0;
end

nNodes = numel(parent.nodes);
nChildren = numel(parent.child);

% indices of original nodes in parent element
iPON=1:3;
% indices of enrichment nodes in parent element
% iPEN=4:nNodes;
Xel=[nodeCoords(parent.nodes(iPON),1),...
     nodeCoords(parent.nodes(iPON),2)]';       
% calculate shape function derivatives of original element
[~, DN] = shape_funct([0;0], 1);
Bel = DN/(Xel*DN);


for c=1:numel(parent.child)
    if (parent.child(c).isTriangle)
        shape = 1;
    else
        shape = 2;
    end
    
    Xch = [nodeCoords(parent.child(c).nodes,1),...
           nodeCoords(parent.child(c).nodes,2)]';

    if (parent.child(c).isTriangle)
        for i = 1:size(gauss.elem,2)
             [N,B,detJ] ...
                = compute_child_element_NBJ_polyIGFEM(iPON,...
                                                      parent.child(c).locPaEnNodes,...
                                                      parent.child(c).locEnNodes, ...
                                                      nNodes, ...
                                                      Xel, ...
                                                      Xch, ...
                                                      parent.child(c).isTriangle, ...
                                                      Bel,...
                                                      gauss.elem(1:end-1,i));                                                  
            factor = detJ*gauss.elem(end,i);
            Nch = shape_funct(gauss.elem(1:end-1,i),shape);
            Xglo = Xch*Nch;
            u_du = refSoln(Xglo(1),Xglo(2));
            UU = UUR(parent.nodes)'*N;
            varargout{1} = varargout{1} + (UU-u_du(1))^2*factor;
            varargout{2} = varargout{2} + u_du(1)^2*factor;
            if (nargout > 2)
                dU = UUR(parent.nodes)'*B;
                varargout{3} = varargout{3} + ((UU - u_du(1))^2 ...
                                            + (dU(1) - u_du(2))^2 ...
                                            + (dU(2) - u_du(3))^2)*factor;
                varargout{4} = varargout{4} + sum(u_du.^2)*factor;
            end
        end
    else
        for i = 1:size(gauss.quadElem,2)
             [N,B,detJ] ...
                = compute_child_element_NBJ_polyIGFEM(iPON,...
                                                      parent.child(c).locPaEnNodes,...
                                                      parent.child(c).locEnNodes, ...
                                                      nNodes, ...
                                                      Xel, ...
                                                      Xch, ...
                                                      parent.child(c).isTriangle, ...
                                                      Bel,...
                                                      gauss.quadElem(1:end-1,i));
            factor = detJ*gauss.quadElem(end,i);
            Nch = shape_funct(gauss.quadElem(1:end-1,i),shape);
            Xglo = Xch*Nch;                                                                      
            u_du = refSoln(Xglo(1),Xglo(2));
            UU = UUR(parent.nodes)'*N;
            varargout{1} = varargout{1} + (UU-u_du(1))^2*factor;
            varargout{2} = varargout{2} + u_du(1)^2*factor;
            if (nargout > 2)
                dU = UUR(parent.nodes)'*B;
                varargout{3} = varargout{3} + ((UU - u_du(1))^2 ...
                                            + (dU(1) - u_du(2))^2 ...
                                            + (dU(2) - u_du(3))^2)*factor;
                varargout{4} = varargout{4} + sum(u_du.^2)*factor;
            end
        end
    end

    
    
 
end






end
