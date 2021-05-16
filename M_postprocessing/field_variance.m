 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Marcus Tan on 11/3/2014
%%% Modified by Marcus Tan
%%% Last modified date: 2/12/2015
%%% Copyright 2014 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if estimate is use, provide updated UUR
function [var,ave,totArea] = field_variance(elem,nodeCoords,UUR,gauss)

nanInd = isnan(UUR);
if (any(nanInd))
	warning('removing NAN values from UUR') 
	UUR(nanInd) = 0.0; % set nan to 0 in case the equation is singular
end

totArea = 0.0;
ave =  0.0;
var = 0.0;
if (isempty(gauss))
    for el = 1:size(elem.elem_node,1)
        if (elem.parent(el).type > 1)
            for c = 1:numel(elem.parent(el).child)
                childArea = abs(polygonArea(nodeCoords(elem.parent(el).child(c).nodes,:)));
                totArea = totArea + childArea;
                val = mean(UUR(elem.parent(el).child(c).nodes));
                ave= ave + val*childArea;
                var = var + val*val*childArea;
            end
        else
            elemArea = abs(polygonArea(nodeCoords(elem.elem_node(el,:),:)));
            totArea = totArea + elemArea;
            val = mean(UUR(elem.elem_node(el,:)));
            ave = ave + val*elemArea;
            var = var + val*val*elemArea;
        end
    end
    ave = ave/totArea;
    var = var/totArea - ave*ave;
    
else
    iPON = 1:3;
    for el = 1:size(elem.elem_node,1)
        Xel = nodeCoords(elem.elem_node(el,:),:)';
        [~,DN] = shape_funct([0;0],1);
        J = Xel*DN;
        Bel = DN/J;
        % polynomial IGFEM
        if (elem.parent(el).type == 2)
            nNodes = numel(elem.parent(el).nodes);
            for c = 1:numel(elem.parent(el).child)
                Xch = nodeCoords(elem.parent(el).child(c).nodes,:)';
                for j = 1:size(gauss.elem,2)
                      [Np2p,~,detJ] ...
                          = compute_child_element_NBJ_polyIGFEM(...
                                iPON, ...
                                elem.parent(el).child(c).locPaEnNodes, ...
                                elem.parent(el).child(c).locEnNodes,...
                                nNodes,...
                                Xel,...
                                Xch,...
                                elem.parent(el).child(c).isTriangle,...
                                Bel,...
                                gauss.elem(1:(end-1),j));
                      val = UUR(elem.parent(el).nodes)'*Np2p;
                      factor = gauss.elem(end,j)*detJ;
                      totArea = totArea + factor;
                      ave = ave + val*factor;
                      var = var + val*val*factor;
                end               
            end
        % NURBS IGFEM   
        elseif (elem.parent(el).type == 3)
            nNodes = numel(elem.parent(el).nodes);
            for c = 1:numel(elem.parent(el).child)
                nurbsSurf = elem.parent(el).child(c).nurbsSurf;
                dnurbsSurf = nrbderiv(nurbsSurf);
                nuknot=length(nurbsSurf.knots{1});
                nvknot=length(nurbsSurf.knots{2});
                uDeg = elem.parent(el).child(c).nurbsSurf.order(1) - 1;
                vDeg = elem.parent(el).child(c).nurbsSurf.order(2) - 1;
                nurbsInd=elem.parent(el).child(c).nurbsInd(elem.parent(el).child(c).locEnNodes);

                nuSpans = nuknot - 2*nurbsSurf.order(1) + 1;
                nvSpans = nvknot - 2*nurbsSurf.order(2) + 1;
                dimensions = [nuSpans,nvSpans,gauss.qua1Dt.npt,gauss.qua1Dn.npt];
                totGaussPts = prod(dimensions);
                rParams = zeros(2,totGaussPts);
                detJp2n = zeros(totGaussPts,1);
                weights = zeros(totGaussPts,1);

                % get all the quantities needed for evaluating the shape functions at
                % all the gauss points
                for i=nurbsSurf.order(1):nuknot-nurbsSurf.order(1)
                   for j=nurbsSurf.order(2):nvknot-nurbsSurf.order(2)
                       rSubEl = [nurbsSurf.knots{1}(i),nurbsSurf.knots{1}(i+1),...
                                 nurbsSurf.knots{1}(i+1),nurbsSurf.knots{1}(i);...
                                 nurbsSurf.knots{2}(j),nurbsSurf.knots{2}(j),...
                                 nurbsSurf.knots{2}(j+1),nurbsSurf.knots{2}(j+1)];
                        for k = 1:gauss.qua1Dt.npt   
                            for m = 1:gauss.qua1Dn.npt
                                linearInd = sub2ind(dimensions,i-uDeg,j-vDeg,k,m);
                                weights(linearInd) = gauss.qua1Dt.weight(k)*gauss.qua1Dn.weight(m);
                                [Np2n,DNp2n]=shape_funct([gauss.qua1Dt.pt(k);gauss.qua1Dn.pt(m)],2);
                                % calculate quantities from parametric space to natural
                                % coordinate space
                                Jp2n=rSubEl*DNp2n;                       
                                detJp2n(linearInd)=det(Jp2n);
                                rParams(:,linearInd)=rSubEl*Np2n;                       
                            end
                        end
                    end
                end
                calcDer = true;
                % compute shape functions
                [Np2p,~,detJp2p] = compute_child_element_NBJ(nurbsSurf, ...
                                                                dnurbsSurf,...
                                                                nurbsInd, ...
                                                                iPON, ...
                                                                elem.parent(el).child(c).locPaEnNodes, ...
                                                                nNodes, ...
                                                                Xel, ...
                                                                Bel, ...
                                                                calcDer, ...
                                                                rParams);  
                factor = detJp2n.*detJp2p.*weights;
                totArea = totArea + sum(factor);
                val = UUR(elem.parent(el).nodes)'*Np2p;
                ave = ave + val*factor;
                var = var + sum(val.*val.*factor');
            
            end
        else
            detJ = det(J);
            for j = 1:size(gauss.elem,2)
                N = shape_funct(gauss.elem(1:(end-1),j), 1);
                factor = gauss.elem(end,j)*detJ;
                totArea = totArea + factor;
                val = UUR(elem.elem_node(el,:))'*N;
                ave = ave + val*factor;
                var = var + val*val*factor;
            end
        end
    end
    ave = ave/totArea;
    var = var/totArea - ave*ave;
end


end