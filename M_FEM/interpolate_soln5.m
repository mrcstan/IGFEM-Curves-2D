%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Marcus Tan on 6/12/2013
%%% Last modified date: 10/27/2014
%%% Copyright 2013 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT:
%   epsge: geometric tolerance for nurbsSurf_phy2param. Default is 1e-8
%   calcDer: logical indicating whether the derivative is calculated.
%            Default is true
% OUTPUT:
%   u_du:  a nx3 array of u, du/dx and du/dy. n is not necessarily the same
%           as length(xI). if a point is not found in any element, then
%           that point will be discarded.
%   xI: a length n  array where u_du is evaluated
%   yI:
%   del: the indices of the points in xI or yI whose nurbs parameters
%        cannot be found.
function [u_du,xI,yI,del]=interpolate_soln(nodeCoords,elem,UUR,xI,yI,epsge,calcDer)
nxI=length(xI);
nyI=length(yI);
if(nxI~=nyI)
    error('interpolate soln: length of xI and yI must be equal')
end
if (nargin < 6)
   epsge = 1e-8; 
end
if (nargin < 7)
    calcDer = true;
end
if (calcDer)
    u_du = NaN(nxI,3);
else
    u_du = NaN(nxI,1);
end

% indices of original nodes in parent element
iPON=1:3;
for j=1:nxI
    foundElem = false;
    for i=1:elem.n_elem
		xElem=nodeCoords(elem.elem_node(i,:),1);
		yElem=nodeCoords(elem.elem_node(i,:),2);   
		IN=inpolygon(xI(j),yI(j),xElem,yElem);
        if(IN)      
            if(elem.parent(i).type == 3)
				% indices of enrichment nodes in parent element
				%iPEN=4:numel(elem.parent(i).nodes);
                for k=1:numel(elem.parent(i).child)
                    %if(inpolygon(xI(j),yI(j),elem.parent(i).child(k).hullx,...
                    %                         elem.parent(i).child(k).hully))  % disabling this check seems faster if the mex function nurbsSurf_phy2param is used                                     
                        rParam = nurbsSurf_phy2param([xI(j);yI(j)],...
                                                     elem.parent(i).child(k).nurbsSurf.knots{1},...
                                                     elem.parent(i).child(k).nurbsSurf.knots{2},...
                                                     elem.parent(i).child(k).nurbsSurf.coefs([1,2,4],:,:),...
                                                     2,...
                                                     epsge);                    
                        if(isempty(rParam))
                            continue
                        end
                        %[~,iCENP,iCEN]=intersect(elem.parent(i).nodes(iPEN),...
                        %                elem.parent(i).child(k).nodes,'stable');
                        %iCENP=iPEN(iCENP);
                        nurbsInd=elem.parent(i).child(k).nurbsInd(elem.parent(i).child(k).locEnNodes);
                        
                                    
                        [Np2p,Bp2p] ...
                                = compute_child_element_NBJ(elem.parent(i).child(k).nurbsSurf, ...
                                                            elem.parent(i).child(k).dnurbsSurf,...
                                                            nurbsInd, ...
                                                            iPON, ...
                                                            elem.parent(i).child(k).locPaEnNodes, ...
                                                            numel(elem.parent(i).nodes),...
                                                            [xElem,yElem]', ...
                                                            [], ...
                                                            calcDer, ...
                                                            rParam);             
                        if (calcDer)    
                            u_du(j,1)=UUR(elem.parent(i).nodes)'*Np2p;
                            u_du(j,2:3)=UUR(elem.parent(i).nodes)'*Bp2p;
                        else
                            u_du(j,1)=UUR(elem.parent(i).nodes)'*Np2p;
                        end
                        foundElem = true;
                        break;
                    %end
                end

            elseif (elem.parent(i).type == 2)
                 for k = 1:numel(elem.parent(i).child)
                    xCh = nodeCoords(elem.parent(i).child(k).nodes,1);
                    yCh = nodeCoords(elem.parent(i).child(k).nodes,2);
                    %[in,~,Xloc] = inTriangle(xI(k),yI(k),xCh,yCh,epsge); 
                    in=inpolygon(xI(j),yI(j),xCh,yCh);
                    if (in)
                        Xch = [xCh,yCh]';
                        if (elem.parent(i).child(k).isTriangle)
                            shape = 1;
                        else
                            shape = 2;
                        end
                        Xloc = local_coord([xI(j),yI(j)]',Xch,shape);
                        [Np2p,Bp2p] ...
                                = compute_child_element_NBJ_polyIGFEM(...
                                        iPON, ...
                                        elem.parent(i).child(k).locPaEnNodes, ...
                                        elem.parent(i).child(k).locEnNodes,...
                                        numel(elem.parent(i).nodes),...
                                        [xElem,yElem]',...
                                        Xch,...
                                        elem.parent(i).child(k).isTriangle,...
                                        [],...
                                        Xloc);
                        u_du(j,1) = UUR(elem.parent(i).nodes)'*Np2p;
                        u_du(j,2:3) = UUR(elem.parent(i).nodes)'*Bp2p;
                        foundElem = true;
                        break
                    end
                end 
            else
                X_local = local_coord([xI(j);yI(j)],...
                                     [xElem';yElem'],1);
                [N, DN] = shape_funct(X_local, 1);
                u_du(j,1)= UUR(elem.elem_node(i,:))'*N;

                if (calcDer)
                    J=[xElem';yElem']*DN;
                    B=DN/J;
                    u_du(j,2:3)=UUR(elem.elem_node(i,:))'*B;
                end
                break;
            end
            
        end
        if (foundElem)
            break
        end
    end
end
% remove NaN's
del = isnan(u_du(:,1));
xI(del) = [];
yI(del) = [];
u_du(del,:) = [];

end