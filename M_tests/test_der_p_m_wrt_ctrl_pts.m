clear all
close all
path(path,'../M_channels')
path(path,'../ChannelFiles')
path(path, '../../NURBS/nurbs_toolbox')
%channelFile = 'DK_DP_check_branchingA_2params_w_diams.channel';
%channelFile = 'DK_DP_check_branchingA_2params.channel';
%channelFile = 'parallel4_start.channel';
channelFile = 'parallel2_start_w_diams.channel';
[channels,designParams] ...
    = read_channels(channelFile);
[designParams,channels.designParamNum] ...
        = design_params2channel_params_map(designParams,...
                                           channels);
xo = zeros(designParams.nParams,1);    
restrictedParams.iniVals = [];

%channels.massin = [];
%channels.powerXdensity = 5.0;
%channels.powerXdensity = [];
[channels,pressureOld,massOld] ...
                        = update_channels(xo,...
                                          designParams, ...
                                          restrictedParams, ...
                                          channels, ...
                                          'add',true);  
figure
plot_channel_network(channels.contvty,channels.nurbs,pressureOld/1000.0,massOld*1e6*60/1065)                                      
D2mass = channels.D2mcf/channels.heatCapacity;
D2pressure = channels.D2P;

del = 1e-7;
method = 'central';
FD_D2pressure = nan(numel(pressureOld),designParams.nParams);
FD_D2mass = nan(numel(massOld),designParams.nParams);

for i = 1:designParams.nParams
    if (strcmpi(method,'forward')) 
        delParams = zeros(designParams.nParams,1);
        delParams(i) = delParams(i) + del; 
        [~,pressureNew1,massNew1] ...
                        = update_channels(delParams,...
                                          designParams, ...
                                          restrictedParams, ...
                                          channels, ...
                                          'add',false);
        FD_D2pressure(:,i) = (pressureNew1 - pressureOld)/del;
        FD_D2mass(:,i) = (massNew1 - massOld)/del;
    elseif (strcmpi(method,'central'))
        delParams = zeros(designParams.nParams,1);
        delParams(i) = delParams(i) + del;  
        [~,pressureNew2,massNew2] ...
                        = update_channels(delParams,...
                                          designParams, ...
                                          restrictedParams, ...
                                          channels, ...
                                          'add',false);
        
        delParams = zeros(designParams.nParams,1);
        delParams(i) = delParams(i) - del; 
        [~,pressureNew1,massNew1] ...
                        = update_channels(delParams,...
                                          designParams, ...
                                          restrictedParams, ...
                                          channels, ...
                                          'add',false);        
        FD_D2pressure(:,i) = 0.5*(pressureNew2 - pressureNew1)/del;
        FD_D2mass(:,i) = 0.5*(massNew2 - massNew1)/del;
    else
        error('unrecognized method')
    end
end

absD2pressureDiff = abs(FD_D2pressure - D2pressure);
relD2pressureDiff = abs(FD_D2pressure./D2pressure - 1);
excludedInd = isnan(relD2pressureDiff) | isinf(relD2pressureDiff); 
absD2massDiff = abs(FD_D2mass - D2mass);
fprintf('pressure grad abs diff = %g\n',max(absD2pressureDiff(:)));
fprintf('pressure grad rel diff = %g\n',max(relD2pressureDiff(~excludedInd)));
fprintf('mass grad abs diff = %g\n',max(absD2massDiff(:)));
                                            