
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Marcus Tan on 8/26/2013
%%% Last modified date: 7/15/2013
%%% Copyright 2013 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this function computes the stiffness matrix of an element split by a
% nurbs interface.
function [KEL,PEL,errflag]=compute_parent_element_nurbs(nodeCoords,...
                          parent,elemHeatSource,heatSourceFunc,convect,itrface,gauss)     

KEL = zeros(parent.nNode, parent.nNode); % Initilize element stifness matrix to zero
PEL = zeros(parent.nNode, 1); % Initilize element force matrix to zero
errflag = false(1,parent.nChild);

% indices of original nodes in parent element
iPON=1:3;
% indices of enrichment nodes in parent element
iPEN=4:parent.nNode;
XOrg=[nodeCoords(parent.node(iPON),1),...
           nodeCoords(parent.node(iPON),2)]';
% for constraint equation
row = zeros(parent.nConstraint,parent.nNode);  
cstrApplied = false(parent.nConstraint,1);
elemHeatSource = elemHeatSource+convect.coef*convect.Tref;

for c=1:parent.nChild

    % iCEN refers to the indices of the enrichment nodes of the child
    % element 
    [~,iCENP,iCEN]=intersect(parent.node(iPEN),parent.child(c).node,'stable');
    iCENP=iPEN(iCENP);
    % iPEN are the indices of the enrichment nodes in the parent element
    nurbsSurf=parent.child(c).nurbsSurf;
    nuknot=length(nurbsSurf.knots{1});
    nvknot=length(nurbsSurf.knots{2});
    dnurbsSurf=parent.child(c).dnurbsSurf;    
    nurbsInd=parent.child(c).nurbsInd(iCEN);
    % make nurbs shape functions
    %{
    if(~parent.child(c).isTriangle && length(parent.child(c).nurbsSeg)>1 ...
        && length(unique(parent.child(c).lineSource))==3)      
        [~,iSCEN] = ismember(parent.child(c).node(iCEN),...
                            [parent.child(c).lineSource(1),...
                            parent.child(c).lineSource(end)]);
       [nurbsShape,dnurbsShape]=create_nurbs_enrichment_function(nurbsSurf,nurbsInd,iSCEN);                 
    else
    %}
    
    [nurbsShape,dnurbsShape]=create_nurbs_enrichment_function(nurbsSurf,nurbsInd);

    % conductivity
    Cmat=parent.child(c).conductivity*eye(2);
    %[Cmat] = constitutive (COMPUTE,conductivity);
    for i=nurbsSurf.order(1):nuknot-nurbsSurf.order(1)
        for j=nurbsSurf.order(2):nvknot-nurbsSurf.order(2)
             if(nurbsSurf.knots{1}(i)~=nurbsSurf.knots{1}(i+1) ...
               && nurbsSurf.knots{2}(j)~=nurbsSurf.knots{2}(j+1))
                rSubEl=[nurbsSurf.knots{1}(i),nurbsSurf.knots{1}(i+1),...
                        nurbsSurf.knots{1}(i+1),nurbsSurf.knots{1}(i);...
                        nurbsSurf.knots{2}(j),nurbsSurf.knots{2}(j),...
                        nurbsSurf.knots{2}(j+1),nurbsSurf.knots{2}(j+1)];
                for k = 1:gauss.qua1Dt.npt
                    for m = 1:gauss.qua1Dn.npt
                        r = [gauss.qua1Dt.pt(k),gauss.qua1Dn.pt(m)];
                        w = gauss.qua1Dt.weight(k)*gauss.qua1Dn.weight(m);
                        % calculate quantities from parametric space to natural
                        % coordinate space
                        [Np2n,DNp2n]=shape_funct(r,2);
                        Jp2n=rSubEl*DNp2n;
                        detJp2n=det(Jp2n);
                        rParam=rSubEl*Np2n;
                        % compute shape functions
                        [Np2p,Bp2p,detJp2p,rPhy] = compute_child_element_NBJ(...
                                                            nurbsSurf, ...
                                                            dnurbsSurf, ...
                                                            nurbsShape, ...
                                                            dnurbsShape,...
                                                            iPON, ...
                                                            iCENP, ...
                                                            parent.nNode, ...
                                                            XOrg, ...
                                                            rParam); 
                        if(detJp2p<=0)
                            errflag(c) = true; 
                        end
                        factor=detJp2n*detJp2p*w;
                        KEL  = KEL + (Bp2p*Cmat*Bp2p'+convect.coef*(Np2p*Np2p'))*factor;
                        PEL = PEL + Np2p*elemHeatSource*factor;
                        if(~isempty(heatSourceFunc))
                            PEL = PEL+Np2p*heatSourceFunc(rPhy(1),rPhy(2))*factor;
                        end
                    end
                end      
             end
        end
    end
    

    % add contribution of line source to stiffness matrix
    % ASSUMPTION: there is at most two line sources per child element
    for i=1:parent.child(c).nLineSource
        nknot=length(parent.child(c).nurbsSeg(i).knots);
        for j=parent.child(c).nurbsSeg(i).order:(nknot-parent.child(c).nurbsSeg(i).order)
            if(parent.child(c).nurbsSeg(i).knots(j)~=parent.child(c).nurbsSeg(i).knots(j+1))
                xiSubElem=[parent.child(c).nurbsSeg(i).knots(j),parent.child(c).nurbsSeg(i).knots(j+1)];
                for k=1:gauss.qua1Dt.npt
                    xi = gauss.qua1Dt.pt(k);
                    w = gauss.qua1Dt.weight(k);
                    % calculate quantities from parametric space to natural
                    % coordinate space
                    [Np2n,DNp2n] = shape_funct_1D(xi);
                    xiParam = xiSubElem*Np2n;
                    Jp2n = xiSubElem*DNp2n;
                    % calculate quantities from physical space to parametric space 
                    rParam = [0;0];
                    rParam(parent.child(c).uv(i)) = parent.child(c).uvParam(i);
                    rParam(setdiff(1:2,parent.child(c).uv(i))) = xiParam;
                     % compute shape functions
                    [Np2p,Bp2p,~,~] = compute_child_element_NBJ(nurbsSurf, ...
                                                               dnurbsSurf, ...
                                                               nurbsShape, ...
                                                               dnurbsShape,...
                                                               iPON, ...
                                                               iCENP, ...
                                                               parent.nNode, ...
                                                               XOrg, ...
                                                               rParam);                   
                    %dnurbsSeg=parent.child(c).dnurbsSeg(i);
                    [~,dxiParam]=nrbdeval(parent.child(c).nurbsSeg(i),parent.child(c).dnurbsSeg(i),xiParam);
                    % fluid flows in the direction increasing parametric coordinate
                    % this model assumes mean temperature = wall
                    % temperature                 
                    % more accurate model
                    if(itrface.modelType == 2)
                        dist2inletPcoord = parent.child(c).lineSourceNP(i,1:2)*Np2n; 
                        dist2inlet = itrface.length(parent.child(c).seg(i))*dist2inletPcoord;
                        conductVal = conductance(dist2inlet,itrface.kapf(parent.child(c).seg(i)),...
                                                itrface.mcf(parent.child(c).seg(i)));
                        factor = norm(dxiParam(1:2))*Jp2n*w; 
                        KEL = KEL + 0.5*conductVal*(Np2p*Np2p')*factor;
                        PEL = PEL + 0.5*conductVal*itrface.Tin(parent.child(c).seg(i))*Np2p*factor;
                    else
                          KEL=KEL+0.5*itrface.mcf(parent.child(c).seg(i))*Np2p*(Bp2p*dxiParam(1:2))'*Jp2n*w; 
                    end
                    
                end
            end
        end
    end
    
    for i = 1:parent.nConstraint
        if (~cstrApplied(i))
            XGlo = nodeCoords(parent.node(parent.cstrLocNode(i)),:)';
            rParam = nurbsSurf_phy2param(XGlo,...
                                     parent.child(c).nurbsSurf.knots{1},...
                                     parent.child(c).nurbsSurf.knots{2},...
                                     parent.child(c).nurbsSurf.coefs([1,2,4],:,:),...
                                     2, 1e-8);

            if(~isempty(rParam))        
                Np2p = compute_child_element_NBJ(nurbsSurf, ...
                                                 dnurbsSurf, ...
                                                 nurbsShape, ...
                                                 dnurbsShape,...
                                                 iPON, ...
                                                 iCENP, ...
                                                 parent.nNode, ...
                                                 XOrg,...
                                                 rParam);

                row(i,:) = Np2p';
                cstrApplied(i) = true;
            end
        end
    end
end


if(any(~cstrApplied))
    disp('compute_parent_element_nurbs: cannot find nurbs parametric coordinates for constraint nodes')
end

%% replace appropriate row (corresponding to constraint node) of stiffness
%  matrix with the constraint equation
for i = 1:parent.nConstraint    
    XGlo = nodeCoords(parent.node(parent.cstrLocNode(i)),:)';
    XLoc = local_coord(1,XGlo,XOrg);
    [N,~] = shape_funct(XLoc,1); 
    row(i,iPON) = N';  
    KEL(parent.cstrLocNode(i),:) = row(i,:);
    PEL(parent.cstrLocNode(i)) = parent.cstrVal(i);
end

end
