%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Marcus Tan on 4/15/2015
%%% Copyright 2015 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this function read the values of a sample from a Latin hypercube sampling
% file that includes the diameter and convert it to one without the
% diameters
clear all
close all
path(path,'../ChannelFiles')
path(path,'../../NURBS/nurbs_toolbox')
channelFile = 'parallel_stephen_start_w_diams.channel'; % channel file with diameters as design parameters
[channels,designParams] = read_channels(channelFile);
diamInds = strcmpi(designParams.type,'DIAM');
nDiamParams = nnz(diamInds);

inputFile = 'parallelStephen_NE_b.lhs'; 

fileID = fopen(inputFile,'r');
line = fgetl(fileID);
if (strncmpi(line,'nSamples',8))
    splitStr = regexp(line,',','split');
    nSamples = str2double(splitStr{2});
else
    error('number of samples must be specified at the beginning of the file')
end
line = fgetl(fileID);
if (strncmpi(line,'nParams',7))
    splitStr = regexp(line,',','split');
    nParams = str2double(splitStr{2});
else
    error('number of parameters must be specified in the second line')
end

if (nParams > designParams.nParams)
    error('number of parameters in sample file is greater than the design parameters in the channel file')
end

allVals = fscanf(fileID,'%g',inf);
if (numel(allVals) ~= nSamples*nParams)
    error('total number of values is inconsistent with number of samples and number of parameters')
end
fclose(fileID);

splitStr = textscan(inputFile,'%s','delimiter','.');
sampleFileType = splitStr{1}{2};

[fileName,pathName] = uiputfile(['*.',sampleFileType],'Save as');
outputFile = [pathName,fileName,];
%outputFile = ['test.',sampleFileType];
fileID = fopen(outputFile,'w');
fprintf(fileID,'nSamples,%i\n',nSamples);
fprintf(fileID,'nParams,%i\n',nParams+nDiamParams);
for i = 1:nSamples
    i;
    j = (i-1)*nParams+1;
    k = j+nParams-1;
    fprintf(fileID,'%g \n',allVals(j:k));
    randDiams = randomized_bounded_values(designParams.bounds(diamInds,:));
    fprintf(fileID,'%g \n',randDiams);
end
fclose(fileID);
