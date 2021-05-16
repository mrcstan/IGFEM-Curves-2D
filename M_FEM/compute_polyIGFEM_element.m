
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Marcus Tan on 8/26/2013
%%% Last modified date: 10/27/2014
%%% Copyright 2013 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this function computes the stiffness matrix of an element split by a
% nurbs interface.
% SUPG can be applied 
function [KEL,PEL,errflag]=compute_polyIGFEM_element(nodeCoords,...
                                                              parent,...
                                                              elemHeatSource,...
                                                              heatSourceFunc,...
                                                              convect,...
                                                              channels,...
                                                              gauss,...
                                                              supg)    
%timings = zeros(8,1);
jacTol = 1e-13;
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



% for constraint equation
nConstraints = numel(parent.cstrLocNodes);
cstrEqns = zeros(nNodes,nConstraints);  
cstrApplied = false(nConstraints,1);
elemHeatSource = elemHeatSource+convect.coef*convect.Tref;


% don't need SUPG for const heat flux model
%if (channels.model == 2)
%    supg = false;
%end

% conductivity
Cmat=parent.conductivity*eye(2);
for c=1:numel(parent.child)
    if (parent.child(c).isTriangle)
        shape = 1;
    else
        shape = 2;
    end
    nChannels = numel(parent.child(c).channelNum);   
    if (nChannels)
        channelVecs = nan(2,nChannels);
        for i = 1:nChannels
            channelVecs(:,i) = (nodeCoords(parent.child(c).channelNodes{i}(end),:) ...
                               -nodeCoords(parent.child(c).channelNodes{i}(1),:))';
        end
        %channelVecs = (nodeCoords(parent.child(c).nodes(parent.child(c).channelLocNodes(2,:)),:) ...
        %              -nodeCoords(parent.child(c).nodes(parent.child(c).channelLocNodes(1,:)),:))';
        if (supg)      
            %channelUnitVecs = bsxfun(@(x,y) x./y,channelVecs,sqrt(sum(channelVecs.^2,1)));
            error('supg not yet implemented for poly IGFEM')
        end
    end
    
    
    Xch = [nodeCoords(parent.child(c).nodes,1),...
           nodeCoords(parent.child(c).nodes,2)]';
    if (supg && nChannels)
        error('supg not yet implemented for poly IGFEM')
    else
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
                KEL = KEL + B*Cmat*B'*factor;
                PEL = PEL + N*elemHeatSource*factor;  
                if(~isempty(heatSourceFunc))
                    Nch = shape_funct(gauss.elem(1:end-1,i),shape);
                    Xglo = Xch*Nch;
                    PEL = PEL + N*heatSourceFunc(Xglo(1),Xglo(2))*factor;
                end
            end
            if (detJ <= jacTol)
                errflag(ch) = true;
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
                KEL = KEL + B*Cmat*B'*factor;
                PEL = PEL + N*elemHeatSource*factor; 
                if(~isempty(heatSourceFunc))
                    Nch = shape_funct(gauss.quadElem(1:end-1,i),shape);
                    Xglo = Xch*Nch;
                    PEL = PEL + N*heatSourceFunc(Xglo(1),Xglo(2))*factor;
                end
                if (detJ <= jacTol)
                    errflag(ch) = true;
                end
            end
        end
    end
    if parent.child(c).isTriangle
        % auxiliary set for determining edge number
        auxSet = [1,2;2,3;1,3;
                  2,1;3,2;3,1];
        nEdges = 3;
    else
        % auxiliary set for determining edge number
        auxSet = [1,2;2,3;3,4;1,4;
                  2,1;3,2;4,3;4,1];
        nEdges = 4;
    end
    if (channels.model == 1) % mean temperature mod
        if (supg && nChannels)
            error('supg not yet implemented for poly IGFEM')
        else
            for i = 1:nChannels
                [~,edge] = ismember(parent.child(c).channelLocNodes(:,i)',auxSet,'rows');
                edge = rem(edge-1,nEdges)+1;
                for j = 1:size(gauss.line,2)
                    Xglo = nodeCoords(parent.child(c).channelNodes{i}(1),:)' ...
                           + channelVecs(:,i)*gauss.line(1,j);
                    Xloc = local_coord_along_edge(Xglo,Xch,shape,edge);
                    %{
                    if(parent.child(c).isTriangle)
                        if (any(Xloc < 0) || any(Xloc > 1) || sum(Xloc) > 1 || any(isnan(Xloc)))
                            warning('triangle: out of range local coord')
                            edge
                            Xch
                            disp(Xglo)
                            disp(Xloc)
                        end
                    else
                        if (any(Xloc < -1) || any(Xloc > 1) || any(isnan(Xloc)))
                            warning('quad: out of range local coord')
                            disp(Xloc)
                        end
                    end
                    %}
                    [N,B] ...
                        = compute_child_element_NBJ_polyIGFEM(iPON,...
                                                              parent.child(c).locPaEnNodes,...
                                                              parent.child(c).locEnNodes, ...
                                                              nNodes, ...
                                                              Xel, ...
                                                              Xch, ...
                                                              parent.child(c).isTriangle, ...
                                                              Bel,...
                                                              Xloc);
                   
                    KEL = KEL + 0.5*channels.mcf(parent.child(c).channelNum(i)) ...
                               *N*(B*channelVecs(:,i))'*gauss.line(end,j);
                end
            end
        end
        
    elseif (channels.model == 2) % constant heat flux model
        error('channel model 2 not yet implemented')
    end
    % enforce Dirichlet BC with Lagrange multiplier method 
    for i = 1:nConstraints
        if (~cstrApplied(i) && any(parent.child(c).locPaNodes == parent.cstrLocNodes(i)))
            Xglo = nodeCoords(parent.nodes(parent.cstrLocNodes(i)),:)';
            Xloc = local_coord(Xglo,Xch,shape);                                                       
            N = compute_child_element_NBJ_polyIGFEM(iPON,...
                                                      parent.child(c).locPaEnNodes,...
                                                      parent.child(c).locEnNodes, ...
                                                      nNodes, ...
                                                      Xel, ...
                                                      Xch, ...
                                                      parent.child(c).isTriangle, ...
                                                      Bel,...
                                                      Xloc);                          

            cstrEqns(:,i) = N;
            cstrApplied(i) = true;
        end
    end
end




%% replace appropriate row (corresponding to constraint node) of stiffness
%  matrix with the constraint equation
if (nConstraints)
    if(any(~cstrApplied))
        warning('cannot find constraint nodes \n')
    else
        KEL = [KEL,zeros(nNodes,nConstraints);zeros(nConstraints,nNodes),zeros(nConstraints,nConstraints)];
        PEL = [PEL;zeros(nConstraints,1)];
        paddedCstrEqns = [cstrEqns;zeros(nConstraints,nConstraints)];
        for i = 1:nConstraints    
            Xglo = nodeCoords(parent.nodes(parent.cstrLocNodes(i)),:)';
            Xloc = local_coord(Xglo,Xel, 1);
            N = shape_funct(Xloc,1); 
            paddedCstrEqns(iPON,i) = N;
            KEL(:,nNodes+i) = paddedCstrEqns(:,i);
            KEL(nNodes+i,:) = paddedCstrEqns(:,i)';
            PEL(nNodes+i) = parent.cstrVals(i);
        end
    end
end

end
