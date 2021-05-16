
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Marcus Tan on 8/26/2013
%%% Last modified date: 10/27/2014
%%% Copyright 2013 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this function computes the stiffness matrix of an element split by a
% nurbs interface.
% SUPG can be applied 
function [KEL,PEL,errflag,timings]=compute_parent_element_nurbs(nodeCoords,...
                                                        parent,...
                                                        elemHeatSource,...
                                                        heatSourceFunc,...
                                                        convect,...
                                                        channels,...
                                                        gauss,...
                                                        supg)     
timings = zeros(2,1);                                                    
nNodes = numel(parent.nodes);
nChildren = numel(parent.child);
nConstraints = numel(parent.cstrLocNode);
KEL = zeros(nNodes, nNodes); % Initilize element stifness matrix to zero
PEL = zeros(nNodes, 1); % Initilize element force matrix to zero
errflag = false(1,nChildren);

% indices of original nodes in parent element
iPON=1:3;
Xel=[nodeCoords(parent.nodes(iPON),1),...
           nodeCoords(parent.nodes(iPON),2)]';
% calculate shape function derivatives of original element
[~, DN] = shape_funct([0;0], 1);
Bel = DN/(Xel*DN);

% calculate nurbs shape function derivatives
calcDer = true;   

% for constraint equation
row = zeros(nConstraints,nNodes);  
cstrApplied = false(nConstraints,1);
elemHeatSource = elemHeatSource+convect.coef*convect.Tref;


% don't need SUPG for const heat flux model
%if (channels.modelType == 2)
%    supg = false;
%end

complement12 = [2,1];
for c=1:numel(parent.child)
    nurbsSurf=parent.child(c).nurbsSurf;
    dnurbsSurf=parent.child(c).dnurbsSurf; 
    nuknot=length(nurbsSurf.knots{1});
    nvknot=length(nurbsSurf.knots{2});
       
    nurbsInd=parent.child(c).nurbsInd(parent.child(c).locEnNodes);
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
    end
    % conductivity
    Cmat=parent.child(c).conductivity*eye(2);
    tic
    %[Cmat] = constitutive (COMPUTE,conductivity);
    for i=nurbsSurf.order(1):nuknot-nurbsSurf.order(1)
        for j=nurbsSurf.order(2):nvknot-nurbsSurf.order(2)
             %if(nurbsSurf.knots{1}(i)~=nurbsSurf.knots{1}(i+1) ...
             %  && nurbsSurf.knots{2}(j)~=nurbsSurf.knots{2}(j+1))
                rSubEl=[nurbsSurf.knots{1}(i),nurbsSurf.knots{1}(i+1),...
                        nurbsSurf.knots{1}(i+1),nurbsSurf.knots{1}(i);...
                        nurbsSurf.knots{2}(j),nurbsSurf.knots{2}(j),...
                        nurbsSurf.knots{2}(j+1),nurbsSurf.knots{2}(j+1)];
                for k = 1:gauss.qua1Dt.npt
                    for m = 1:gauss.qua1Dn.npt
                        % calculate quantities from parametric space to natural
                        % coordinate space
                        [Np2n,DNp2n]=shape_funct([gauss.qua1Dt.pt(k);gauss.qua1Dn.pt(m)],2);
                        Jp2n=rSubEl*DNp2n;
                        detJp2n=det(Jp2n);
                        rParam=rSubEl*Np2n;
                        % compute shape functions
                        [Np2p,Bp2p,detJp2p,rPhy] ...
                             = compute_child_element_NBJ(nurbsSurf, ...
                                                         dnurbsSurf,...
                                                         nurbsInd, ...
                                                         iPON, ...
                                                         parent.child(c).locPaEnNodes, ...
                                                         nNodes, ...
                                                         Xel, ...
                                                         Bel, ...
                                                         calcDer, ...
                                                         rParam); 
    
                        
                        % SUPG weight function
                        Wfunc = Np2p;
                        %{
                        if (supg)
                            % strategy 1
                            %
                            for n = 1:nChannels
                                [he,Bsw] = streamwise_elem_length(unitTan(n,:)',Bp2p);
                                Wfunc = Wfunc+0.5*he*Bsw;
                            end
                            %
                            % strategies 2 and 3
                            %{
                            [he,Bsw] = streamwise_elem_length(unitTan',Bp2p);
                            Wfunc = Wfunc+0.5*he*Bsw;
                            %}
                        end
                        %}
                        
                        factor=detJp2n*detJp2p*gauss.qua1Dt.weight(k)*gauss.qua1Dn.weight(m);
                        KEL = KEL + (Bp2p*Cmat*Bp2p'+convect.coef*(Wfunc*Np2p'))*factor;
                        PEL = PEL + Wfunc*elemHeatSource*factor;
                        if(~isempty(heatSourceFunc))
                            PEL = PEL+Wfunc*heatSourceFunc(rPhy(1),rPhy(2))*factor;
                        end
                        if(~errflag(c) && detJp2p<=0)
                            errflag(c) = true; 
                        end
                    end
                end      
             %end
        end
    end
    timings(1) = timings(1) + toc;

    % add contribution of line source to stiffness matrix
    % ASSUMPTION: there is at most two line sources per child element
    tic
    for i=1:nChannels
        nknots=length(parent.child(c).nurbsSeg(i).knots);
        for j=parent.child(c).nurbsSeg(i).order:(nknots-parent.child(c).nurbsSeg(i).order)
            %if(parent.child(c).nurbsSeg(i).knots(j)~=parent.child(c).nurbsSeg(i).knots(j+1))
                xiSubElem=[parent.child(c).nurbsSeg(i).knots(j),parent.child(c).nurbsSeg(i).knots(j+1)];
                for k=1:gauss.qua1Dt.npt
                    % calculate quantities from parametric space to natural
                    % coordinate space
                    [Np2n,DNp2n] = shape_funct_1D( gauss.qua1Dt.pt(k));
                    xiParam = xiSubElem*Np2n;
                    Jp2n = xiSubElem*DNp2n;
                    % calculate quantities from physical space to parametric space 
                    rParam([parent.child(c).uv(i), ...
                            complement12(parent.child(c).uv(i))]) ...
                            = [parent.child(c).uvParam(i),xiParam];
                     % compute shape functions

                      [Np2p,Bp2p] ...
                             = compute_child_element_NBJ(nurbsSurf, ...
                                                         dnurbsSurf,...
                                                         nurbsInd, ...
                                                         iPON, ...
                                                         parent.child(c).locPaEnNodes, ...
                                                         nNodes, ...
                                                         Xel, ...
                                                         Bel, ...
                                                         calcDer, ...
                                                         rParam);                                        
                    % SUPG weight function                                       
                    Wfunc = Np2p; 
                    %{
                    if (supg)
                        % strategy 1
                        %
                        for n = 1:nChannels
                            [he,Bsw] = streamwise_elem_length(unitTan(n,:)',Bp2p);
                            Wfunc = Wfunc+0.5*he*Bsw;
                        end
                        %
                        % strategies 2 and 3
                        %{
                        [he,Bsw] = streamwise_elem_length(unitTan',Bp2p);
                        Wfunc = Wfunc+0.5*he*Bsw;
                        %}
                    end
                    %}
                    %dnurbsSeg=parent.child(c).dnurbsSeg(i);
                    [~,dxiParam]=nrbdeval(parent.child(c).nurbsSeg(i),parent.child(c).dnurbsSeg(i),xiParam);
                    % fluid flows in the direction increasing parametric coordinate
                    % this model assumes mean temperature = wall
                    % temperature                 
                    % more accurate model
                    if(channels.modelType == 2)
                        dist2inletPcoord = parent.child(c).channelNurbsParam(1:2,i)'*Np2n; 
                        dist2inlet = channels.length(parent.child(c).channelNum(i))*dist2inletPcoord;
                        conductVal = conductance(dist2inlet,channels.kapf(parent.child(c).channelNum(i)),...
                                                channels.mcf(parent.child(c).channelNum(i)));
                        factor = norm(dxiParam(1:2))*Jp2n*gauss.qua1Dt.weight(k); 
                        KEL = KEL + 0.5*conductVal*(Wfunc*Np2p')*factor;
                        PEL = PEL + 0.5*conductVal*channels.Tin(parent.child(c).channelNum(i))*Wfunc*factor;
                    else
                          KEL=KEL+0.5*channels.mcf(parent.child(c).channelNum(i))...
                                     *Wfunc*(Bp2p*dxiParam(1:2))'...
                                     *Jp2n*gauss.qua1Dt.weight(k); 
                    end
                    
                end
            %end
        end
    end
    timings(2) = timings(2) + toc;
    
    for i = 1:nConstraints
        if (~cstrApplied(i))
            XGlo = nodeCoords(parent.nodes(parent.cstrLocNode(i)),:)';
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
                                                 parent.child(c).locPaEnNodes, ...
                                                 nNodes, ...
                                                 Xel,...
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
for i = 1:nConstraints    
    XGlo = nodeCoords(parent.nodes(parent.cstrLocNode(i)),:)';
    XLoc = local_coord(1,XGlo,Xel);
    [N,~] = shape_funct(XLoc,1); 
    row(i,iPON) = N';  
    KEL(parent.cstrLocNode(i),:) = row(i,:);
    PEL(parent.cstrLocNode(i)) = parent.cstrVal(i);
end

end
