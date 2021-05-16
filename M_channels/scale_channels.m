close all
clear
path(path,'../../NURBS/nurbs_toolbox')
channelFile = '../ChannelFiles/parallel5x5ref.channel';
shift = [0,0];
scaling = 2.5;

splitStr = regexp(channelFile,'/','split');
headers = {['# scaled ',splitStr{end}]};
[channels,designParams] = read_channels(channelFile);

options.showChannels = false;
options.showPts = false;
options.showDesignParams = true;
options.showParamLabels = true;

% scale and shift channel end points and NURBS control points
channels.pts = bsxfun(@plus,channels.pts*scaling,shift); 
for i = 1:numel(channels.nurbs)
    channels.nurbs(i).coefs(1:2,:) = bsxfun(@plus,channels.nurbs(i).coefs(1:2,:)*scaling,shift');
end

% scale and shift design parameter bounds
designParams.bounds = bsxfun(@plus,designParams.bounds*scaling,shift);


plot_channels_and_design_parameters(channels.contvty, ...
                                    channels.pts, ...
                                    channels.nurbs,...
                                    designParams,options)
[filename,pathname] = uiputfile('*.channel','Save as','../ChannelFiles/');
outFile = [pathname,filename];
write_channel_file(outFile,channels,designParams,'w',headers)



