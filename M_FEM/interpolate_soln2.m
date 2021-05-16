%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Marcus Tan on 6/12/2013
%%% Last modified date: 10/15/2013
%%% Copyright 2013 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUT:
%   u_du:  a nx3 array of u, du/dx and du/dy. n is not necessarily the same
%           as length(xI). if a point is not found in any element, then
%           that point will be discarded.
%   xI: a length n  array where u_du is evaluated
%   yI:
%   del: the indices of the points in xI or yI whose nurbs parameters
%        cannot be found.
function [u_du,xI,yI,del]=interpolate_soln(nodeCoord,elem,UUR,xI,yI)
nxI=length(xI);
nyI=length(yI);
if(nxI~=nyI)
    error('interpolate soln: length of xI and yI must be equal')
end
u_du=NaN(nxI,3);
shape=1; 
for j=1:nxI
    for i=1:elem.n_elem
		xElem=nodeCoord(elem.elem_node(i,:),1);
		yElem=nodeCoord(elem.elem_node(i,:),2);   
		IN=inpolygon(xI(j),yI(j),xElem,yElem);
        if(IN)      
            if(elem.parent(i).isParent)
				% indices of original nodes in parent element
				iPON=1:3;
				% indices of enrichment nodes in parent element
				iPEN=4:elem.parent(i).nNode;
                for k=1:elem.parent(i).nChild
                    %if(inpolygon(xI(j),yI(j),elem.parent(i).child(k).hullx,...
                    %                         elem.parent(i).child(k).hully))  % disabling this check seems faster if the mex function nurbsSurf_phy2param is used                                     
                        rParam = nurbsSurf_phy2param([xI(j);yI(j)],...
                                                     elem.parent(i).child(k).nurbsSurf.knots{1},...
                                                     elem.parent(i).child(k).nurbsSurf.knots{2},...
                                                     elem.parent(i).child(k).nurbsSurf.coefs([1,2,4],:,:),...
                                                     2);                    
                        if(isempty(rParam))
                            continue
                        end
                        [~,iCENP,iCEN]=intersect(elem.parent(i).node(iPEN),...
                                        elem.parent(i).child(k).node,'stable');
                        iCENP=iPEN(iCENP);
                        nurbsInd=elem.parent(i).child(k).nurbsInd(iCEN);
                             % make nurbs shape functions
                        if(~elem.parent(i).child(k).isTriangle ...
                            && length(elem.parent(i).child(k).nurbsSeg)>1 ...
                            && length(unique(elem.parent(i).child(k).lineSource))==3)
                            [~,iSCEN] = ismember(elem.parent(i).child(k).node(iCEN),...
                                                [elem.parent(i).child(k).lineSource(1),...
                                                elem.parent(i).child(k).lineSource(end)]);
                            [nurbsShape,dnurbsShape]...
                                    =create_nurbs_enrichment_function(...
                                            elem.parent(i).child(k).nurbsSurf,nurbsInd,iSCEN);                
                        else    
                            [nurbsShape,dnurbsShape]...
                            =create_nurbs_enrichment_function(...
                            elem.parent(i).child(k).nurbsSurf,nurbsInd);
                        end
                        [Np2p,Bp2p,~,~]=compute_child_element_NBJ(...
                                        elem.parent(i).child(k).nurbsSurf,...
                                        elem.parent(i).child(k).dnurbsSurf,...
                                        nurbsShape,dnurbsShape,...
                                        iPON,iCENP,elem.parent(i).nNode,[xElem,yElem]',rParam); 
                        u_du(j,1)=UUR(elem.parent(i).node)'*Np2p;
                        u_du(j,2:3)=UUR(elem.parent(i).node)'*Bp2p;
                        break;
                    %end
                end
			else
				X_local = local_coord(shape, [xI(j),yI(j)],...
									 [xElem';yElem']);
				[N, DN] = shape_funct(X_local, shape);
				J=[xElem';yElem']*DN;
				B=DN/J;
				u_du(j,1)= UUR(elem.elem_node(i,:))'*N;
				u_du(j,2:3)=UUR(elem.elem_node(i,:))'*B;
				break;
            end
        end
    end
end
% remove NaN's
del = isnan(u_du(:,1));
xI(del) = [];
yI(del) = [];
u_du(del,:) = [];

end