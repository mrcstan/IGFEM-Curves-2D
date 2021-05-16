
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Marcus Tan on 8/26/2013
%%% Last modified date: 10/27/2014
%%% Copyright 2013 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this function computes the stiffness matrix of an element split by a
% nurbs interface.
% SUPG can be applied 
function [KEL,PEL,errflag]=compute_parent_element_nurbs(nodeCoords,...
                                                              parent,...
                                                              elemHeatSource,...
                                                              heatSourceFunc,...
                                                              convect,...
                                                              channels,...
                                                              gauss,...
                                                              supg)    
%timings = zeros(8,1);

nNodes = numel(parent.nodes);
nChildren = numel(parent.child);
KEL = zeros(nNodes, nNodes); % Initialize element stifness matrix to zero
PEL = zeros(nNodes, 1); % Initialize element force matrix to zero
errflag = false(1,nChildren);

% indices of original nodes in parent element
iPON=1:3;
% indices of enrichment nodes in parent element
% iPEN=4:nNodes;
Xel=[nodeCoords(parent.nodes(iPON),1),...
           nodeCoords(parent.nodes(iPON),2)]';       
% calculate shape function derivatives of original element
[~, DN] = shape_funct([0;0], 1);
Bel = DN/(Xel*DN);

% calculate nurbs shape function derivatives
calcDer = true;

% for constraint equation
nConstraints = numel(parent.cstrLocNodes);
cstrEqns = zeros(nNodes,nConstraints);  
cstrApplied = false(nConstraints,1);
elemHeatSource = elemHeatSource+convect.coef*convect.Tref;


% don't need SUPG for const heat flux model
%if (channels.model == 2)
%    supg = false;
%end

complement12 = [2,1];
for c=1:numel(parent.child)
    nChannels = numel(parent.child(c).channelNum);    
    if (supg && nChannels)
        unitTan = nodeCoords(parent.child(c).nodes(parent.child(c).channelLocNodes(2,:)),:)...
                  -nodeCoords(parent.child(c).nodes(parent.child(c).channelLocNodes(1,:)),:);      
        % strategies 2 and 3    
        %{
        % unitTan = bsxfun(@rdivide,unitTan,sqrt(sum(unitTan.^2,2))); % stragegy 3
        if (parent.nLineSource > 1)
            unitTan = sum(bsxfun(@times,unitTan,channels.mcf(parent.child(c).seg)),1);
        end
        %}
        unitTan = unitTan/norm(unitTan,2); % IS THIS RIGHT????????????????
        %unitTan = bsxfun(@(x,y) x./y,unitTan,sqrt(sum(unitTan.^2,2)));
    end
   
    
    % conductivity
    Cmat=parent.conductivity*eye(2);
    %[Cmat] = constitutive (COMPUTE,conductivity);
    
    nurbsSurf = parent.child(c).nurbsSurf;
    dnurbsSurf = nrbderiv(nurbsSurf);
    nuknot=length(nurbsSurf.knots{1});
    nvknot=length(nurbsSurf.knots{2});
    uDeg = parent.child(c).nurbsSurf.order(1) - 1;
    vDeg = parent.child(c).nurbsSurf.order(2) - 1;
    nurbsInd=parent.child(c).nurbsInd(parent.child(c).locEnNodes);
    
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
    %rPhy = nrbeval(nurbsSurf,{rParams(1,:),rParams(2,:)});
                                               
    % compute shape functions
    [Np2p,Bp2p,detJp2p,rPhy] = compute_child_element_NBJ(nurbsSurf, ...
                                                         dnurbsSurf,...
                                                         nurbsInd, ...
                                                         iPON, ...
                                                         parent.child(c).locPaEnNodes, ...
                                                         nNodes, ...
                                                         Xel, ...
                                                         Bel, ...
                                                         calcDer, ...
                                                         rParams);                                                     
    if(any(detJp2p<=0))
        errflag(c) = true; 
    end
    factor=detJp2n.*detJp2p.*weights;

    if (supg)
        for i = 1:totGaussPts
            j = 2*i-1;
            % SUPG weight function
            Wfunc = Np2p(:,i);
            % strategy 1
            for n = 1:nChannels
                [he,Bsw] = streamwise_elem_length(unitTan(n,:)',Bp2p(:,j:j+1));
                Wfunc = Wfunc+0.5*he*Bsw;
            end
            %
            % strategies 2 and 3
            %{
            [he,Bsw] = streamwise_elem_length(unitTan',Bp2p);
            Wfunc = Wfunc+0.5*he*Bsw;
            %}
            KEL = KEL + (Bp2p(:,j:j+1)*Cmat*Bp2p(:,j:j+1)'+convect.coef*(Wfunc*Np2p(:,i)'))*factor(i);
            PEL = PEL + Wfunc*elemHeatSource*factor(i);
            if(~isempty(heatSourceFunc))
                PEL = PEL + Wfunc*heatSourceFunc(rPhy(1,i),rPhy(2,i))*factor(i);
            end
        end
    else
        for i = 1:totGaussPts
            j = 2*i-1;
            KEL = KEL + (Bp2p(:,j:j+1)*Cmat*Bp2p(:,j:j+1)'+convect.coef*(Np2p(:,i)*Np2p(:,i)'))*factor(i);
            PEL = PEL + Np2p(:,i)*elemHeatSource*factor(i);
            if(~isempty(heatSourceFunc))
                PEL = PEL + Np2p(:,i)*heatSourceFunc(rPhy(1,i),rPhy(2,i))*factor(i);
            end
        end
    end
    
   

    % add contribution of line source to stiffness matrix
    % ASSUMPTION: there is at most two line sources per child element
    %
    for i=1:nChannels
        nknots = length(parent.child(c).nurbsSeg(i).knots);
        deg = parent.child(c).nurbsSeg(i).order - 1;
        nSpans = nknots - 2*parent.child(c).nurbsSeg(i).order + 1;
        dimensions = [nSpans,gauss.qua1Dt.npt];
        totGaussPts = prod(dimensions);
        rParams = zeros(2,totGaussPts);
        xiParams = zeros(totGaussPts,1);
        Jp2n = zeros(totGaussPts,1);
        weights = zeros(totGaussPts,1);
        dist2inlet = zeros(totGaussPts,1); % for constant heat flux model
       
        for j=parent.child(c).nurbsSeg(i).order:(nknots-parent.child(c).nurbsSeg(i).order)
            xiSubElem=[parent.child(c).nurbsSeg(i).knots(j),parent.child(c).nurbsSeg(i).knots(j+1)];
            for k=1:gauss.qua1Dt.npt
                linearInd = sub2ind(dimensions,j-deg,k);
                weights(linearInd) = gauss.qua1Dt.weight(k);
                % calculate quantities from parametric space to natural
                % coordinate space
                [Np2n,DNp2n] = shape_funct_1D(gauss.qua1Dt.pt(k));
                xiParams(linearInd) = xiSubElem*Np2n;
                Jp2n(linearInd) = xiSubElem*DNp2n;
                % calculate quantities from physical space to parametric space 
                rParams(parent.child(c).uv(i),linearInd) = parent.child(c).uvParam(i);
                rParams(complement12(parent.child(c).uv(i)),linearInd) = xiParams(linearInd); 
                dist2inlet(linearInd) = channels.length(parent.child(c).channelNum(i))...
                                        *parent.child(c).channelNurbsParam(1:2,i)'*Np2n;
            end
        end
        %timings(4) = timings(4) + toc;
         % compute shape functions
        [Np2p,Bp2p] = compute_child_element_NBJ(nurbsSurf, ...
                                                dnurbsSurf,...
                                                nurbsInd, ...
                                                iPON, ...
                                                parent.child(c).locPaEnNodes, ...
                                                nNodes, ...
                                                Xel, ...
                                                Bel, ...
                                                calcDer, ...
                                                rParams);

        %dnurbsSeg=parent.child(c).dnurbsSeg(i);
        [~,dxiParams]=nrbdeval(parent.child(c).nurbsSeg(i),parent.child(c).dnurbsSeg(i),xiParams);
        % fluid flows in the direction increasing parametric coordinate
        if(channels.model == 2)
            for j = 1:totGaussPts
                conductVal = conductance(dist2inlet(j),channels.kapf(parent.child(c).channelNum(i)),...
                                         channels.mcf(parent.child(c).channelNum(i)));
                factor = norm(dxiParams(1:2,j))*Jp2n(j)*weights(j); 
                KEL = KEL + 0.5*conductVal*(Np2p(:,j)*Np2p(:,j)')*factor(j);
                PEL = PEL + 0.5*conductVal*channels.Tin(parent.child(c).channelNum(i))*Np2p(:,j)*factor(j);
            end
        else
            if (supg)
                for j = 1:totGaussPts
                    k = 2*j-1;
                    Wfunc = Np2p(:,j);
                    % strategy 1
                    for n = 1:nChannels
                        [he,Bsw] = streamwise_elem_length(unitTan(n,:)',Bp2p(:,k:k+1));
                        Wfunc = Wfunc+0.5*he*Bsw;
                    end
                    %
                    % strategies 2 and 3
                    %{
                    [he,Bsw] = streamwise_elem_length(unitTan',Bp2p);
                    Wfunc = Wfunc+0.5*he*Bsw;
                    %}

                    KEL=KEL+0.5*channels.mcf(parent.child(c).channelNum(i))...
                               *Wfunc*(Bp2p(:,k:k+1)*dxiParams(1:2,j))'...
                               *Jp2n(j)*weights(j); 
                end
            else
                for j = 1:totGaussPts
                    k = 2*j-1;
                    KEL=KEL+0.5*channels.mcf(parent.child(c).channelNum(i))...
                               *Np2p(:,j)*(Bp2p(:,k:k+1)*dxiParams(1:2,j))'...
                               *Jp2n(j)*weights(j); 
                end
            end
        end
        %timings(6) = timings(6) + toc;
    end
    %
    % enforce Dirichlet BC with Lagrange multiplier method 
    for i = 1:nConstraints
        if (~cstrApplied(i) && any(parent.child(c).locPaNodes == parent.cstrLocNodes(i)))
            Xglo = nodeCoords(parent.nodes(parent.cstrLocNodes(i)),:)';
            rParam = nurbsSurf_phy2param(Xglo,...
                                         parent.child(c).nurbsSurf.knots{1},...
                                         parent.child(c).nurbsSurf.knots{2},...
                                         parent.child(c).nurbsSurf.coefs([1,2,4],:,:),...
                                         2, 1e-8);

            if(~isempty(rParam))                                                     
                Np2p = compute_child_element_NBJ(nurbsSurf, ...
                                                 dnurbsSurf,...
                                                 nurbsInd, ...
                                                 iPON, ...
                                                 parent.child(c).locPaEnNodes, ...
                                                 nNodes, ...
                                                 Xel, ...
                                                 Bel, ...
                                                 false, ...
                                                 rParam);                             

                cstrEqns(:,i) = Np2p;
                cstrApplied(i) = true;
            end
        end
    end
end




%% replace appropriate row (corresponding to constraint node) of stiffness
%  matrix with the constraint equation
if (nConstraints)
    if(any(~cstrApplied))
        fprintf('cannot find nurbs parametric coordinates for constraint nodes \n')
    else
        KEL = [KEL,zeros(nNodes,nConstraints);zeros(nConstraints,nNodes),zeros(nConstraints,nConstraints)];
        PEL = [PEL;zeros(nConstraints,1)];
        paddedCstrEqns = [cstrEqns;zeros(nConstraints,nConstraints)];
        for i = 1:nConstraints    
            Xglo = nodeCoords(parent.nodes(parent.cstrLocNodes(i)),:)';
            Xloc = local_coord(Xglo,Xel,1);
            N = shape_funct(Xloc,1); 
            paddedCstrEqns(iPON,i) = N;
            KEL(:,nNodes+i) = paddedCstrEqns(:,i);
            KEL(nNodes+i,:) = paddedCstrEqns(:,i)';
            PEL(nNodes+i) = parent.cstrVals(i);
        end
    end
end

end
