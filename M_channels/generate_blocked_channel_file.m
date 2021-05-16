close all
clear
path(path,'../channelFiles')
path(path,'../../NURBS/nurbs_toolbox')
path(path, '../M_geom_toolbox')
path(path, '../M_preFEM')
path(path,'../../MatlabUsefulFunctions/export_fig')
channelFile = '3x3ref_deg6.channel';

nBlocksPerSet = 1; % max number of blocked channels per set
[channels,designParams] = read_channels(channelFile);

[filename,pathname] = uiputfile('*.blk','Save as','../blockedChannelFiles/');
outFile = [pathname,filename];
fileID = fopen(outFile,'w');

interiorChannels = find(channels.contvty(:,1) ~= channels.inletEndPoint ...
                        & channels.contvty(:,2) ~= channels.pressureOutletEndPoint);
blockedSets = combnk(interiorChannels,nBlocksPerSet);     
formatStr = [repmat('%i\t',1,nBlocksPerSet),'\n'];
fprintf(fileID,formatStr,blockedSets');
fclose(fileID);
