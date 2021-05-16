close all
clear
path(path,'../../NURBS/nurbs_toolbox')
sampleFile = '../SampleFiles/parallel5x5_b.rand';
channelFile = '../ChannelFiles/parallel5x5ref.channel';

shift = [0,0];
scaling = 2.5;

[~,designParams] = read_channels(channelFile);

[~,nSamples,nParams] = read_sample_file(sampleFile,1);

splitStr = strsplit(sampleFile,'.');
[filename,pathname] = uiputfile(['*.',splitStr{end}],'Save as','../SampleFiles/');
outFile = [pathname,filename];

fileID = fopen(outFile,'w');
fprintf(fileID,'nSamples,%i\n',nSamples);
fprintf(fileID,'nParams,%i\n',nParams);

if nParams ~= designParams.nParams
    error('sample file does not have the same number of design parameters as the channel file')
end
for i = 1:nSamples
    iniVals = read_sample_file(sampleFile,i); 
    for j = 1:numel(iniVals)
        if strcmpi(designParams.type{j},'CTRL_PT')
            newVal = iniVals(j)*scaling+shift(designParams.ctrlPtDim(j));
            fprintf(fileID,'%g \n',newVal);
        end
    end
end
fclose(fileID);


               



