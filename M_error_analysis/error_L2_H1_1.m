%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Marcus Tan on 7/26/2013
%%% Last modified date: 8/30/2013
%%% Copyright 2013 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUT:  
% totalL2: totalL2 error
% totalH1: totalH1 error
% totL2norm: totalL2 norm
% totH1norm: totalH1 norm
% elem.parent.L2: for all elements, including regular elements, the
% normalized L2 error 
% elem.parent.H1:
% elem.parent.child.L2: for parent elements, the normalized error of the 
% child elements  
% elem.parent.child.H1:
function [totalL2,totalH1,totParentL2,totParentH1,totL2norm,totH1norm,elem]...
        =error_L2_H1(analytical,UUR,nodeCoord,elem,gauss)
normalized = false;
totalL2 = 0;
totalH1 = 0;
totParentL2 = 0;
totParentH1 = 0;
totL2norm = 0.0;
totH1norm = 0.0;
for elem_num=1:elem.n_elem
        elem.parent(elem_num).L2=0;
        elem.parent(elem_num).H1=0;
        sumNormL2=0;
        sumNormH1=0;
        if(elem.parent(elem_num).isParent)

            for c=1:elem.parent(elem_num).nChild
              [elem.parent(elem_num).child(c).L2,normL2,...
               elem.parent(elem_num).child(c).H1,normH1]...
                  =error_squared_L2_H1_child_element(analytical,UUR,...
                   nodeCoord,elem.parent(elem_num),gauss,c);
               elem.parent(elem_num).L2=elem.parent(elem_num).L2...
                                +elem.parent(elem_num).child(c).L2;
               sumNormL2=sumNormL2+normL2;
               if(normalized)
                    elem.parent(elem_num).child(c).L2 = sqrt(elem.parent(elem_num).child(c).L2/normL2);
               else
                    elem.parent(elem_num).child(c).L2 = sqrt(elem.parent(elem_num).child(c).L2);
               end
               elem.parent(elem_num).child(c).normL2=sqrt(normL2);
               elem.parent(elem_num).H1=elem.parent(elem_num).H1+elem.parent(elem_num).child(c).H1;
               sumNormH1=sumNormH1+normH1;
               if(normalized)
                    elem.parent(elem_num).child(c).H1=sqrt(elem.parent(elem_num).child(c).H1/normH1);
               else
                    elem.parent(elem_num).child(c).H1=sqrt(elem.parent(elem_num).child(c).H1);
               end
               elem.parent(elem_num).child(c).normH1=sqrt(normH1);
            end
            totParentL2 = totParentL2+elem.parent(elem_num).L2;
            totParentH1 = totParentH1+elem.parent(elem_num).H1;
        else
            [elem.parent(elem_num).L2,sumNormL2,....
             elem.parent(elem_num).H1,sumNormH1]...
                  =error_squared_L2_H1_regular_element(analytical,UUR,...
                  nodeCoord,elem.elem_node(elem_num,:),gauss);
         
        end
        totalL2=totalL2+elem.parent(elem_num).L2;
        totalH1=totalH1+elem.parent(elem_num).H1;
        if(normalized)
            elem.parent(elem_num).L2=sqrt(elem.parent(elem_num).L2/sumNormL2);
        else
            elem.parent(elem_num).L2=sqrt(elem.parent(elem_num).L2);
        end
        elem.parent(elem_num).normL2=sqrt(sumNormL2);
        if(normalized)
            elem.parent(elem_num).H1=sqrt(elem.parent(elem_num).H1/sumNormH1);
        else
            elem.parent(elem_num).H1=sqrt(elem.parent(elem_num).H1);
        end
        elem.parent(elem_num).normH1=sqrt(sumNormH1);
        totL2norm = totL2norm + sumNormL2;
        totH1norm = totH1norm + sumNormH1;
end

totalL2=sqrt(totalL2);
totalH1=sqrt(totalH1);
totL2norm = sqrt(totL2norm);
totH1norm = sqrt(totH1norm);
totParentL2 = sqrt(totParentL2);
totParentH1 = sqrt(totParentH1);

end