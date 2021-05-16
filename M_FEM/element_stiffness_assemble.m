%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Ahmad Raeisi Najafi on 10/08/2012
%%% Modified by Marcus Tan
%%% Last modified date: 8/4/2014
%%% Copyright 2012 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This function assembles the element stiffness matrices into global
% matricies, partitioned according to the equation numbering.
% INPUT:
% elem_node: the current element node conectivity
function [Kpp,Kpf,Kfp,Kff,Pp,Pf] = element_stiffness_assemble ...
        (elem_node,eq_num,Kel,Pel,Kpp,Kpf,Kfp,Kff,Pp,Pf)

gloInds = eq_num(elem_node);
pLocInds = gloInds < 0;
fLocInds = gloInds > 0;
pGloInds = -gloInds(pLocInds);
fGloInds = gloInds(fLocInds);

%{
gloInds = combvec(pGloInds',pGloInds');
Kpp =  spsubsasgn(Kpp,gloInds(1,:),gloInds(2,:),Kel(pLocInds,pLocInds),@plus);

gloInds = combvec(pGloInds',fGloInds');
Kpf =  spsubsasgn(Kpf,gloInds(1,:),gloInds(2,:),Kel(pLocInds,fLocInds),@plus);

gloInds = combvec(fGloInds',pGloInds');
Kfp =  spsubsasgn(Kfp,gloInds(1,:),gloInds(2,:),Kel(fLocInds,pLocInds),@plus);

gloInds = combvec(fGloInds',fGloInds');
Kff =  spsubsasgn(Kff,gloInds(1,:),gloInds(2,:),Kel(fLocInds,fLocInds),@plus);
%}

%
Kpp(pGloInds,pGloInds) =  Kpp(pGloInds,pGloInds) + Kel(pLocInds,pLocInds);
Kpf(pGloInds,fGloInds) =  Kpf(pGloInds,fGloInds) + Kel(pLocInds,fLocInds);
Kfp(fGloInds,pGloInds) =  Kfp(fGloInds,pGloInds) + Kel(fLocInds,pLocInds);
Kff(fGloInds,fGloInds) =  Kff(fGloInds,fGloInds) + Kel(fLocInds,fLocInds);
%
Pp(pGloInds) = Pp(pGloInds) + Pel(pLocInds);
Pf(fGloInds) = Pf(fGloInds) + Pel(fLocInds);

end % end of function