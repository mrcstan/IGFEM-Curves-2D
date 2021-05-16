% Note: Tested on MATLAB 2015b
%       graph and maxflow functions are not available in MATLAB 2014a or older
clear
close all
path(path,'../ChannelFiles')
channelFile = 'parallel2x2ref.channel';
channels = read_channels(channelFile);
nChannels = size(channels.contvty,1);
%nNodes = max(channels.contvty(:));
%G = sparse(channels.contvty(:,1),channels.contvty(:,2),ones(nChannels,1),nNodes,nNodes);
G = graph(channels.contvty(:,1),channels.contvty(:,2),ones(nChannels,1));
[maxFlow,flowMatrix] = maxflow(G,channels.inletEndPoint+1, ...
                                      channels.pressureOutletEndPoint-1)
[maxFlow,flowMatrix] = maxflow(G,channels.pressureOutletEndPoint-1,channels.inletEndPoint+1)