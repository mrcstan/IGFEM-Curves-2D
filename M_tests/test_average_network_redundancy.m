clear
close all
path(path,'../ChannelFiles')
path(path,'../M_channels')
path(path,'../../NURBS/nurbs_toolbox')
channelFile = 'parallel_stephen_start.channel';
channels = read_channels(channelFile);
ave = average_network_redundancy(channels,...
                                [channels.inletEndPoint; ...
                                 channels.pressureOutletEndPoint]);
fprintf('average network redundancy = %g \n',ave)