%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Marcus Tan in 2014
%%% Last modified on 10/14/2015
%%% Copyright 2014 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% updates the following variables of the channels
% i) connectivity point coordinates
% ii) control points
% iii) diameters
% iv) lengths
% v)  volume
% vi) branching points 
% vii) kinks 
% viii) if pressure or mass output is requested, also output mass flow rate
%       times heat capacities
% action: 'add': add delParams to initial values of the design parameters
%         'replace': replace initial values of the design parameters with
%                     delParams
%         'updateMassPressure': just update the mass and pressure without
%                               updating the geometry of the channels
% calcGrad: request gradients of the pressure and mass wrt to lengths and
%           diameters be calculated
function  [channels,pressure,mass] = update_channels(delParams,...
                                                     designParams, ...
                                                     restrictedParams, ...
                                                     iniChannels, ...
                                                     action, ...
                                                     calcGrad)


switch action
    case 'add'
         delParams = delParams + [designParams.iniVals;restrictedParams.iniVals]; % new design parameter value
         updateGeom = true;
         updateMassPressure = true;
    case 'replace'
         updateGeom = true;
         updateMassPressure = true;
    case 'updateMassPressure'
         updateGeom = false;
         updateMassPressure = true;         
    otherwise
         error('unknown action')
end

channels = iniChannels;

if updateGeom
    for i = 1:designParams.nParams;
        if (strcmpi(designParams.type{i},'CTRL_PT'))
            dim = designParams.ctrlPtDim(i);
            for j = 1:numel(designParams.channelNum{i})            
                chanNum = designParams.channelNum{i}(j);
                cptInd = designParams.ctrlPtNum{i}(j);
                if dim == 4
                    channels.nurbs(chanNum).coefs(dim,cptInd) = delParams(i);
                else
                    channels.nurbs(chanNum).coefs(dim,cptInd)...
                        = delParams(i)*channels.nurbs(chanNum).coefs(4,cptInd);
                end
            end
        elseif(strcmpi(designParams.type{i},'DIAM'))
            channels.diams(designParams.channelNum{i}) = delParams(i); 
        else
            error('unknown design parameter type')
        end
    end

    for i = 1: channels.nNurbs
        channels.pts(channels.contvty(i,1),:) = channels.nurbs(i).coefs(1:2,1)' ...
                                                 /channels.nurbs(i).coefs(4,1);
        channels.pts(channels.contvty(i,2),:) = channels.nurbs(i).coefs(1:2,end)' ...
                                                 /channels.nurbs(i).coefs(4,end);
    end

    %Compute the length of each channels
    ngpts = 10;
    [channels.length,channels.lineSegs,channels.segLengths]...
                = nurbs_channel_lengths(channels,ngpts);   

    channels.branch_kinks = branching_points_n_kinks(channels);

    if (isfield(channels,'vertexCoords') && ~isempty(channels.vertexCoords) ...
        && isfield(designParams,'vertices2params'))
        for i = 1:size(channels.vertexCoords,1)
            ind = ~isnan(designParams.vertices2params(i,:));     
            channels.vertexCoords(i,ind) = delParams(designParams.vertices2params(i,ind));       
        end
    end
end

if updateGeom || updateMassPressure 
     channels.vol = channel_volume(channels.length,...
                                   channels.diams,...
                                   channels.heights,...
                                   channels.crossSection);
end

if nargout > 1 || updateMassPressure 
    if calcGrad
        [pressure,mass,DLpressure,DLmass,DdiamPressure,DdiamMass] ...
            = network_pressure_mass_flow_rate(channels.contvty,...
                                              channels.nurbs,...
                                              channels.diams,...
                                              channels.heights,...
                                              channels.viscosity,...
                                              channels.inletEndPoint,...
                                              channels.massin,...
                                              channels.powerXdensity,...
                                              channels.pressureOutletEndPoint,...
                                              channels.pressureOut,...
                                              channels.crossSection);
        [D2L, D2diam] = D2_channel_lengths_n_diameters(channels,designParams);
        channels.D2mcf = channels.heatCapacity*(DLmass*D2L + DdiamMass*D2diam);
        channels.D2P = DLpressure*D2L + DdiamPressure*D2diam; % does not include the effect of change in viscosity with temperature                                  
    else
        [pressure,mass] ...
            = network_pressure_mass_flow_rate(channels.contvty,...
                                              channels.nurbs,...
                                              channels.diams,...
                                              channels.heights,...
                                              channels.viscosity,...
                                              channels.inletEndPoint,...
                                              channels.massin,...
                                              channels.powerXdensity,...
                                              channels.pressureOutletEndPoint,...
                                              channels.pressureOut,...
                                              channels.crossSection); 
       channels.D2mcf = zeros(numel(mass),designParams.nParams);
       channels.D2P = zeros(numel(pressure),designParams.nParams);
    end
    channels.mcf = mass*channels.heatCapacity;
  
end
end