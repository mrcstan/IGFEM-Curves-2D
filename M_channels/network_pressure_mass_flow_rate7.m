%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Marcus Tan on 7/10/2014
%%% Last modified date: 10/13/2015
%%% Copyright 2014 University of Illinois at Urbana-Champaign. 
%%% All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this function calculates the pressure at each nodal point of a network of
% channels and the mass flow rate in each channel described by a geom.
% The mass flow rate is positive if the flow is from the starting node to
% the end node.
% FORMULA:
%       mass flow rate = pi*diams^4/(128*kinematic_viscosity*length)...
%                        *pressure drop
%       volumetric flow rate = pi*diams^4/(128*kinematic_viscosity*length)
%                               *g*head_loss
%       pressure drop = density*g*head_loss
% NOTE: the factors for square is 1/28.46 compared to pi/128 for sphere
% INPUT:
%   chan_nodes: a number of channels x 2 matrix of the channels and their
%               end nodes
%   geom: an array containing the geom description of each channel. it can be: 
%         i) an array of nurbs structures or
%         ii) number of nodes x number of dimensions matrix of nodal
%             coordinates
%         iii) an array of the lengths of the channels 
%   diams: an array of the diameter of each channel. If only one value is
%         specified, all channels have the same diameter equal to the
%         single value. For 'circular' and 'square' cross section, this is
%         the diameter and length of a side, respectively. For 'rectangle'
%         cross section, this is the width of the cross section, i.e., the
%         width in the plane of the panel
%   heights: only relevant for rectangle cross section, it is the thickness
%           perpendicular to the plane of the panel           
%   nu: kinematic viscosity of the fluid
%   nodalSources: an array that specifies the nodes where the source/sink
%                  strength is to be specified.
%   sourceStrengths: an array the same length of nodalSources that specifies
%                   the source/sink strength at the corresponding nodes
%                   indicated in nodalSources. A source is given a positive
%                   value.    
%   powerXdensity: power (energy/per unit time)*density (mass/per volume)
%                  For this to take effect, sourceStrengths should be
%                  empty. Also, only nodalSources should only contain one
%                  entry
%   nodalBCs: an array that specifies the nodes where the pressure is to be
%             specified
%   BCs: an array the same length as nodalBCs that specifies
%                  the pressures at nodes specified in nodalBCs 
%   crossSection: 'circular','square','rectangle'
% OUTPUT:
%   pressure: pressure at every node
%   chanMass: the mass flow rate in every channel. the flow rate is
%             positive if the flow is from the starting node to the ending nodes 
%   DLpressure: derivative of the pressure wrt length of each channel
%   DLchanMass: derivative of the mass flow rate wrt length of each channel
%   DdiamPressure: derivative of the pressure wrt diameter of each channel
%   DdiamMass: derivative of the mass flow rate wrt diameter of each
%               channel
%   DnuPressure: derivative of the pressure wrt viscosity
%   DnuMass: derivative of the mass flow rate wrt viscosity
function [pressure, chanMass,...
          DLpressure,DLchanMass,...
          DdiamPressure,DdiamMass,...
          DnuPressure] ...
                = network_pressure_mass_flow_rate(chan_nodes,...
                                                  geom,...
                                                  diams,...
                                                  heights,...
                                                  nu,...
                                                  nodalSources,...
                                                  sourceStrengths,...
                                                  powerXdensity,...
                                                  nodalBCs,...
                                                  BCs,...
                                                  crossSection)
                                              
nChans = size(chan_nodes,1);
nNodes = max(chan_nodes(:));

if numel(diams) == 1
    diams = diams*ones(nChans,1);
elseif numel(diams) ~= nChans
    error(['A single diameter value or diameter value corresponding to each',...
           ' channel should be specified'])
end

if numel(nodalSources) == numel(sourceStrengths) && numel(sourceStrengths) >= 1
    isStandard = true;
elseif isempty(sourceStrengths) && numel(nodalSources) == 1 ...
        && numel(powerXdensity)
    isStandard = false;
else
    error(['if sourceStrengths is nonempty, it must have same number of entries as nodalSources',...
          'Otherwise, one entry for both nodalSources and powerXdensity should be specified'])
end

if numel(nodalBCs) ~= numel(BCs)
     error(['number of values in BCs should correspond to number',...
            ' of nodes in nodalBCs'])
end

if numel(nodalBCs) == 0
    error('At least one nodal BC must be specified') 
end

if any(intersect(nodalSources,nodalBCs))
    error('cannot specify a node both as a source and a BC node')
end

diamTol = (max(diams)-min(diams))*1e-8;
% calculate conductivities of each channel
if (strcmpi(crossSection,'circular'))
    factors = pi/(128.0*nu)*diams.^4; % circular cross section
elseif (strcmpi(crossSection,'square'))
    factors = 1/(28.264*nu)*diams.^4; % square cross section
elseif (strcmpi(crossSection,'rectangular'))
    eps = heights./diams;
    factors = (0.25*(1/3.0 - 64*eps/pi^5.*tanh(0.5*pi./eps))/nu).*heights.^3.*diams;
    zeroInds = (abs(diams) < diamTol) | (abs(heights) < diamTol);
    factors(zeroInds) = 0.0;
else
    error('unknown cross section')
end
ngpts = 10;
conductivities = zeros(nChans,1);
if (isstruct(geom))
    if (numel(geom) ~= nChans)
        error('number of nurbs should be equal to the number of rows in chan_nodes')
    end
    if (isfield(geom,'length'))
        for i = 1:nChans
            conductivities(i) = factors(i)/geom.length(i);
        end
    else       
        len = zeros(nChans,1);
        for i = 1:nChans
            len(i) = nurbs_arc_length(ngpts,geom(i));
            conductivities(i) = factors(i)/len(i);
        end
    end
elseif (isnumeric(geom))
    if (size(geom,2) == 2 || size(geom,2) == 3)
        if (size(geom,1) ~= nNodes)
            error('number of rows of geom should equal to the number of nodes')
        end
        len = sqrt(sum((geom(chan_nodes(:,2),:)-geom(chan_nodes(:,1),:)).^2,2));
        conductivities = factors./len;
    elseif (size(geom,2) == 1)
        if (size(geom,1) ~= nChans)
            error('number of rows of geom should equal number of rows in chan_nodes')
        end
        len = geom;
        conductivities = factors./geom;       
    end
else
     error('unknown geom variable type')
end


% setup matrix
rows = chan_nodes(:,[1,2]);
cols = chan_nodes(:,[2,1]);
diags = [chan_nodes(:,1),chan_nodes(:,2)];
K = sparse([rows,diags],...
           [cols,diags],...
           [-conductivities;-conductivities;conductivities;conductivities],...
           nNodes,nNodes);
%Kb4BC = full(K)*9.81
% impose BC
K(nodalBCs,:) = 0.0; % zero the rows specified in nodalBCs
linearInd = sub2ind([nNodes,nNodes],nodalBCs,nodalBCs);
[rows,cols] = meshgrid(1:nNodes,nodalBCs);
Krhs = sparse(rows(:),...
              cols(:),...
              full(-K(:,nodalBCs))',...
              nNodes,nNodes);
%tempKrhs = full(Krhs)*9.81
Krhs(linearInd) = 1.0;
rhstemp = zeros(nNodes,1);
rhstemp(nodalBCs) = BCs;


K(:,nodalBCs) = 0.0; % zero the columns specified in nodalBCs
K(linearInd) = 1.0; % set the diagonal element corresponding to nodalBCs 

if isStandard
    % set nodal source strenghts
    rhs = zeros(nNodes,1);
    rhs(nodalSources) = sourceStrengths;
    rhs = rhs + Krhs*rhstemp; % transfer nodal pressure BCs to the right hand side   
    pressure = K\rhs;
else
    if (numel(nodalSources) > 1 || numel(nodalBCs) > 1)
        error('specifying power x density only implemented for one nodal source and one nodal BC')
    end
    rhsBC = Krhs*rhstemp;
    K = full(K); % faster to deal with full matrix for small system size
    func = @(x) funcPressure(x,K,rhsBC,nodalSources,powerXdensity,BCs);
    x0 = zeros(nNodes,1);
    x0(nodalSources)  = BCs + 1; % ensures that denominator does not vanish
    x0(nodalBCs) = BCs;
    
    fsolveOptions = optimoptions('fsolve','Algorithm','trust-region-dogleg',...
                                          'DerivativeCheck','off',...
                                          'Jacobian','on',...
                                          'Display','off',...
                                          'TolFun',1e-15);
     % Note: the gradient calculation uses the fact that K is symmetric
     [pressure,~,exitFlag] = fsolve(func,x0,fsolveOptions);
     if exitFlag < 1
         if exitFlag < 1
            error('pressure: exit flag = %i',exitFlag)
        end
     end
end

chanMass = conductivities.*(pressure(chan_nodes(:,1))-pressure(chan_nodes(:,2)));

if (nargout > 2)    
    % derivatives wrt lengths
    % K*P = Min
    % K*DP + DK*P = DMin = 0 (0 because inlet flow rate is fixed)
    DLKProws = zeros(2*nChans,1);
    DLKPcols = [1:nChans;1:nChans];
    DLKPvals = zeros(2*nChans,1);
    for i = 1:nChans
        k1 = 2*i - 1;
        k2 = 2*i;
        DLKProws(k1:k2) = chan_nodes(i,:)'; 
        val = -conductivities(i)/len(i)*(pressure(chan_nodes(i,1))-pressure(chan_nodes(i,2)));
        DLKPvals(k1:k2) = [val;-val];
    end
    DLKP = sparse(DLKProws,DLKPcols(:),DLKPvals,nNodes,nChans);
    % Do not need to set the row of DLKP corresponding to nodalBCs to zero
    % as the unknowns are decoupled from the rest of the other equations 
    % (the diagonal entries in K are set to one to achieve this). 
    % The solutions corresponding to these nodalBCs can then be set to zero
    if ~isStandard
        K(nodalSources,nodalSources) ...
            = K(nodalSources,nodalSources) ...
              + powerXdensity/(pressure(nodalSources)-pressure(nodalBCs))^2; 
    end
    
    DLpressure = full(K\(-DLKP)); % [DP1/DL1,DP1/DL2,...;
                                          %  DP2/DL1,DP2/DL2,...;
                                          %  ....               ]   
    DLpressure(nodalBCs,:) = 0.0; % derivative of prescribed nodal pressure is 0
    DLchanMass = zeros(nChans,nChans); % [Dm1/DL1,Dm1/DL2,...;
                                       %  Dm2/DL1,Dm2/DL2,...;
                                       %  ....               ]
    % mass flow rate = conductivity*pressure drop
    for i = 1:nChans
        DLchanMass(:,i) = conductivities.*(DLpressure(chan_nodes(:,1),i) ...
                                          -DLpressure(chan_nodes(:,2),i));
        DLchanMass(i,i) = DLchanMass(i,i) - conductivities(i)/len(i) ...
                                            *(pressure(chan_nodes(i,1)) ...
                                             - pressure(chan_nodes(i,2)));
    end
    
    % derivative of conductivities wrt to diameters
    if (strcmpi(crossSection,'circular') ...
            || strcmpi(crossSection,'square'))
        DdiamConductivities = 4*conductivities./diams;
    elseif (strcmpi(crossSection,'rectangular'))
        eps = heights./diams;
        DdiamConductivities = 0.25/nu*(1/3.0 - 32/pi^4.*sech(0.5*pi./eps)).*heights.^3./len;
        DdiamConductivities(zeroInds) = 0.0;
    end
    % derivatives wrt diameter
    DdiamKProws = zeros(2*nChans,1);
    DdiamKPcols = [1:nChans;1:nChans];
    DdiamKPvals = zeros(2*nChans,1);
  
    for i = 1:nChans
        k1 = 2*i - 1;
        k2 = 2*i;
        DdiamKProws(k1:k2) = chan_nodes(i,:)'; 
        if diams(i) > diamTol
            val = DdiamConductivities(i)*(pressure(chan_nodes(i,1))-pressure(chan_nodes(i,2)));
        else
            val = 0.0;
        end
        DdiamKPvals(k1:k2) = [val;-val];
    end
    DdiamKP = sparse(DdiamKProws,DdiamKPcols(:),DdiamKPvals,nNodes,nChans);
    DdiamPressure = full(K\(-DdiamKP)); % [DP1/Ddiam1,DP1/Ddiam2,...;
                                            %  DP2/Ddiam1,DP2/Ddiam2,...;
                                            %  ....               ]   
    DdiamPressure(nodalBCs,:) = 0.0; % derivative of prescribed nodal pressure is 0
    DdiamMass = zeros(nChans,nChans); % [Dm1/Ddiam1,Dm1/Ddiam2,...;
                                       % Dm2/Ddiam1,Dm2/Ddiam2,...;
                                       %  ....               ]
    for i = 1:nChans
        DdiamMass(:,i) = conductivities.*(DdiamPressure(chan_nodes(:,1),i) ...
                                          -DdiamPressure(chan_nodes(:,2),i));
        if diams(i) > diamTol
            DdiamMass(i,i) = DdiamMass(i,i) + DdiamConductivities(i) ...
                                                *(pressure(chan_nodes(i,1)) ...
                                                 - pressure(chan_nodes(i,2)));
        end
    end
    % derivative wrt to kinematic viscosity
    DnuKP = zeros(nNodes,1);
    for i = 1:nChans
        val = -conductivities(i)/nu*(pressure(chan_nodes(i,1))-pressure(chan_nodes(i,2)));
        DnuKP(chan_nodes(i,1)) = DnuKP(chan_nodes(i,1)) + val;
        DnuKP(chan_nodes(i,2)) = DnuKP(chan_nodes(i,2)) - val;
    end
    DnuPressure =  full(K\(-DnuKP));
    DnuPressure(nodalBCs) = 0.0;
  
end
end

function [vals,grad] = funcPressure(x,K,rhsBC,nodalSource,powerXdensity,BCs)
    vals =  K*x - rhsBC;
    vals(nodalSource) = vals(nodalSource) ...
                            - powerXdensity/(x(nodalSource)-BCs);
    if nargout > 1
        grad = K; % it should be the transpose but since K is symmetric, K=K'
        grad(nodalSource,nodalSource) = grad(nodalSource,nodalSource) ...
                                       + powerXdensity/(x(nodalSource)-BCs)^2;
                                   
    else
        grad = [];
    end
end

%{
function [vals,grad] = funcDdPressure(x,pressures,K,DdKP, ...
                                      nodalSources,powerXdensity,nodalBCs)
    vals = K*x + DdKP;
    vals(nodalSources) = vals(nodalSources) ...
                        + powerXdensity/(pressures(nodalSources)-pressures(nodalBCs))^2 ...
                         *x(nodalSources);
    if nargout > 1
        grad = K;
        grad(nodalSource,nodalSource) = grad(nodalSource,nodalSource) ...
                            + powerXdensity/(pressures(nodalSources)-pressures(nodalBCs))^2;
    else
        grad = [];
    end
end
%}