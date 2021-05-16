clear all
close all
path(path,'../ChannelFiles')
channelFile = 'parallel2x2ref.channel';
channels = read_channels(channelFile);
nChannels = size(channels.contvty,1);
nNodes = max(channels.contvty(:));
G = sparse(channels.contvty(:,1),channels.contvty(:,2),ones(nChannels,1),nNodes,nNodes);
[maxFlow,flowMatrix] = graphmaxflow(G,channels.inletEndPoint+1, ...
                                      channels.pressureOutletEndPoint-1)
[maxFlow,flowMatrix] = graphmaxflow(G,channels.pressureOutletEndPoint-1,channels.inletEndPoint+1)