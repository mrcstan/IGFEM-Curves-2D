%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Marcus Tan in 2014
%%% Last modified on 2/2/2015
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
function  [channels,pressure,mass] = update_channels(delParams,...
                                                     designParams, ...
                                                     restrictedParams, ...
                                                     iniChannels, ...
                                                     action)

if (strcmpi(action,'add'))
    delParams = delParams + [designParams.iniVals;restrictedParams.iniVals]; % new design parameter value
elseif (strcmpi(action,'replace'))
    
else
    error('unknown action')
end



channels = iniChannels;

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


%nSegs = 100;
%channels.curvePts = cell(channels.nNurbs,1);
for i = 1: channels.nNurbs
    %channels.box_corners{i} = nurbs_curve_box_corners(channels.nurbs(i));
    %[~,channels.spanCurva{i}] = nurbs_curve_total_curvature(ngpts,channels.nurbs(i));
    %channels.spanCurva{i} = channels.spanCurva{i}(:,3);
    %channels.curvePts{i} = linearize_nurbs(channels.nurbs(i),nSegs);
    channels.pts(channels.contvty(i,1),:) = channels.nurbs(i).coefs(1:2,1)' ...
                                             /channels.nurbs(i).coefs(4,1);
    channels.pts(channels.contvty(i,2),:) = channels.nurbs(i).coefs(1:2,end)' ...
                                             /channels.nurbs(i).coefs(4,end);
end
    
%Compute the length of each channels
ngpts = 10;
%channels.length = nurbs_channel_lengths(channels,ngpts);
[channels.length,channels.lineSegs,channels.segLengths]...
            = nurbs_channel_lengths(channels,ngpts);   
channels.vol = channel_volume(channels.length,channels.diams,channels.crossSection);

channels.junc = branching_points(channels.contvty);   
channels.kinks = kinks(channels.nurbs);

if (isfield(channels,'vertexCoords') && ~isempty(channels.vertexCoords) ...
    && isfield(designParams,'vertices2params'))
    for i = 1:size(channels.vertexCoords,1)
        ind = ~isnan(designParams.vertices2params(i,:));     
        channels.vertexCoords(i,ind) = delParams(designParams.vertices2params(i,ind));       
    end
end

if (nargout > 1)
    [pressure,mass,DLpressure,DLmass,DdiamPressure,DdiamMass] ...
        = network_pressure_mass_flow_rate(channels.contvty,...
                                          channels.nurbs,...
                                          channels.diams,...
                                          channels.viscosity,...
                                          channels.inletEndPoint,...
                                          channels.massin,...
                                          channels.pressureOutletEndPoint,...
                                          channels.pressureOut,...
                                          channels.crossSection);                                  
    channels.mcf = mass*channels.heatCapacity;
    [D2L, D2diam] = D2_channel_lengths_n_diameters(channels,designParams);
    %channels.D2mcf = channels.heatCapacity*(DLmass*D2L);
    %channels.D2P = DLpressure*D2L;
    channels.D2mcf = channels.heatCapacity*(DLmass*D2L + DdiamMass*D2diam);
    channels.D2P = DLpressure*D2L + DdiamPressure*D2diam; % does not include the effect of change in viscosity with temperature
end
end