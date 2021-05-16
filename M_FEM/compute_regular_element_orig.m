%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Marcus Tan on 6/11/2013
%%% Last modified date: 7/24/2014
%%% Copyright 2013 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% computer stiffness matrix of regular element where the line source along
% the edge of the element can be described by NURBS curve.
% Allows SUPG to be applied to mean temperature model.
% SUPG cannot be applied to const heat flux model.
function [KEL,PEL, errflag]   =   compute_regular_element(nodeCoords,...
                                                          elem_node,...
                                                          parent,...
                                                          elemHeatSource,...
                                                          conductivity,...
                                                          heatSourceFunc,...
                                                          convect,...
                                                          channels,...
                                                          gauss,...
                                                          supg)
% initalize quantities of interest
KEL = zeros(3, 3);
PEL = zeros(3, 1);
errflag = false;
Xel = [nodeCoords(elem_node,1), ...
       nodeCoords(elem_node,2)]';

shape = 1;

[~, DN] = shape_funct([0;0], shape);
J = Xel*DN; % Xel is a 2x3 matrix, DN is a 3x2 matrix
detJ = J(1,1)*J(2,2)-J(1,2)*J(2,1);
if(detJ<=0)
    warning('Jacobian is non-positive')
    errflag = true;
end
B = DN/J; 

if (channels.modelType == 2)
    supg = false;
end

Bsw = zeros(3,1);
% 3 strategies for applying SUPG to multiple line sources

if (parent.type == 1)
    nChannels = numel(parent.channelNum);
    channelVecs = nodeCoords(parent.channelNodes(2,:),:)-nodeCoords(parent.channelNodes(1,:),:);
    if (supg)
        % Strategy 1
        for i = 1:nChannels            
            [he,Bsw1] = streamwise_elem_length(channelVecs(i,:)'/norm(channelVecs(i,:),2),B);
            Bsw = Bsw+0.5*he*Bsw1;
        end
        %
        % Strategies 2 and 3
        %{
        %unitTan = bsxfun(@rdivide,lineSourceVec,sqrt(sum(lineSourceVec.^2,2)));
        unitTan = lineSourceVec; % strategy 2
        if (parent.nLineSource > 1)
            unitTan = sum(bsxfun(@times,unitTan,channels.mcf(parent.lineSourceSeg)),1);

        end
        unitTan = unitTan/norm(unitTan,2);
        [he,Bsw] = streamwise_elem_length(unitTan',B);
        Bsw = 0.5*he*Bsw;
        %}
    end
end



Cmat=conductivity*eye(2);
%[Cmat] = constitutive (COMPUTE, conductivity);
% add constant part of convective heat transfer to distributed heat source
elemHeatSource = elemHeatSource+convect.coef*convect.Tref;

for j = 1:gauss.tri2D.npt
    r = gauss.tri2D.pt(1:2,j);
    w = gauss.tri2D.weight(j);    
    [N, ~] = shape_funct(r', shape);
    Wfunc = N+Bsw;
    KEL  = KEL + (B*Cmat*B'+convect.coef*(Wfunc*N'))*detJ*w;
    PEL  = PEL + Wfunc*elemHeatSource*detJ*w;
    XX=Xel*N;
    if(~isempty(heatSourceFunc))
        PEL=PEL + Wfunc*heatSourceFunc(XX(1),XX(2))*detJ*w;
    end
end

%
% add contribution of line source to local stiffness matrix
if(parent.type == 1)
    if(channels.modelType == 2)
        for i=1:nChannels
            if(lparent.channelLocEdge(i) == 1)
                dS=sqrt(J(1,1)^2+J(2,1)^2);
            elseif(parent.channelLocEdge(i) == 2)
                dS=sqrt((J(1,2)-J(1,1))^2+(J(2,2)-J(2,1))^2);
            elseif(parent.channelLocEdge(i) ==3)
                dS=sqrt(J(1,2)^2+J(2,2)^2);
            end;
            
            r = zeros(2,1);
            for j = 1:gauss.tri1D.npt
                if(locEdge==1)
                    r=[gauss.tri1D.pt(j);0];
                elseif(locEdge==2)
                    r=[gauss.tri1D.pt(j);1-gauss.tri1D.pt(j)];
                elseif(locEdge==3)
                    r=[0;gauss.tri1D.pt(j)];
                end
                [N, ~] = shape_funct(r, shape);

                %
                % a channel line source
                % since two elements share a line source at the edge, the
                % contribution to KEL is halved.
                % this model assumes mean temperature = wall temperature

                [N1D,~ ] = shape_funct_1D(gauss.tri1D.pt(j));
                dist2inletPcoord = parent.channelNurbsParam(1:2,i)*N1D;
                dist2inlet = channels.length(parent.channelNum(i))*dist2inletPcoord;
                conductVal = conductance(dist2inlet,channels.kapf(parent.channelNum(i)),...
                                         channels.mcf(parent.channelNum(i)));
                factor = dS*gauss.tri1D.weight(j);        
                KEL = KEL + 0.5*conductVal*(N*N')*factor;
                PEL = PEL + 0.5*conductVal*channels.Tin(parent.channelNum(i))*N*factor;

            end
        end
    else
        for i=1:nChannels
            r = zeros(2,1);
            for j = 1:gauss.tri1D.npt
                if(parent.channelLocEdge(i) == 1)
                    r=[gauss.tri1D.pt(j);0];
                elseif(parent.channelLocEdge(i) == 2)
                    r=[gauss.tri1D.pt(j);1-gauss.tri1D.pt(j)];
                elseif(parent.channelLocEdge(i) == 3)
                    r=[0;gauss.tri1D.pt(j)];
                end
                N = shape_funct(r, shape);

                %
                % a channel line source
                % since two elements share a line source at the edge, the
                % contribution to KEL is halved.
                KEL = KEL+0.5*channels.mcf(parent.channelNum(i))*(N+Bsw)...
                             *channelVecs(i,:)*B'*gauss.tri1D.weight(j);
            end
        end
    end
end
end

