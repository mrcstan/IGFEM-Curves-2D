%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Marcus Tan on 6/12/2013
%%% Last modified date: 1/26/2016
%%% Copyright 2013 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT:
%   epsge: geometric tolerance for nurbsSurf_phy2param. Default is 1e-8
%   calcDer: logical indicating whether the derivative is calculated.
%            Default is true
% OUTPUT:
%   u:  a nx3 array of u, where n is not necessarily the same
%           as length(xI). if a point is not found in any element, then
%           that point will be discarded.
%   u_x: a nx3 array of du/dx 
%   u_y: a nx3 array of du/dy
%   xI: a length n  array where u, u_x, u_y is evaluated
%   yI:
%   del: the indices of the points in xI or yI whose nurbs parameters
%        cannot be found.
function [u,u_x,u_y,xI,yI,del]=interpolate_soln(nodeCoords,elem,UUR,xI,yI,epsge,calcDer)
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
u = nan(nxI,1);
if (calcDer)
    u_x = nan(nxI,1);
    u_y = nan(nxI,1);
else
    u_x = [];
    u_y = [];
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
                            u(j)=UUR(elem.parent(i).nodes)'*Np2p;
                            u_x(j)=UUR(elem.parent(i).nodes)'*Bp2p(:,1);
                            u_y(j)=UUR(elem.parent(i).nodes)'*Bp2p(:,2);
                        else
                            u(j)=UUR(elem.parent(i).nodes)'*Np2p;
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
                        u(j) = UUR(elem.parent(i).nodes)'*Np2p;
                        u_x(j) = UUR(elem.parent(i).nodes)'*Bp2p(:,1);
                        u_y(j) = UUR(elem.parent(i).nodes)'*Bp2p(:,2);
                        foundElem = true;
                        break
                    end
                end 
            else
                X_local = local_coord([xI(j);yI(j)],...
                                     [xElem';yElem'],1);
                [N, DN] = shape_funct(X_local, 1);
                u(j)= UUR(elem.elem_node(i,:))'*N;

                if (calcDer)
                    J=[xElem';yElem']*DN;
                    B=DN/J;
                    u_x(j)=UUR(elem.elem_node(i,:))'*B(:,1);
                    u_y(j)=UUR(elem.elem_node(i,:))'*B(:,2);
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
del = isnan(u);
xI(del) = [];
yI(del) = [];
u(del) = [];
if calcDer
    u_x(del) = [];
    u_y(del) = [];
end
end