%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Marcus Tan on 2/5/2015
%%% Modified on 2/15/2015
%%% Copyright 2015 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this function read the values of an initial design from an initial design
% file
% FILE FORMAT:
%{
nSamples, <number of initial designs>
nParams, <number of design parameters>
val11
val12
...
val1nParams
...
...
valnSamples1
valnSamples2
...
valnSamplesnParams
%}
% FILE EXAMPLES:
%   See the files with extensions lhs,rand,rand1,smmp in the ChannelFiles
%   directory
function [vals,nSamples,nParams] = read_sample_file(inputFile,sampleNum)
fileID = fopen(inputFile,'r');
line = fgetl(fileID);
splitStr = regexp(line,',','split');
splitStr{1} = strtrim(splitStr{1});
if (strcmpi(splitStr{1},'nSamples'))
    nSamples = str2double(splitStr{2});
    if (sampleNum > nSamples)
         error('sample number not found in sample file')
    end
else
    error('number of samples must be specified at the beginning of the file')
end
line = fgetl(fileID);
splitStr = regexp(line,',','split');
splitStr{1} = strtrim(splitStr{1});
if (strcmpi(splitStr{1},'nParams'))
    nParams = str2double(splitStr{2});
else
    error('number of parameters must be specified in the second line')
end

allVals = fscanf(fileID,'%g',inf);
if (numel(allVals) ~= nSamples*nParams)
    error('total number of values is inconsistent with number of samples and number of parameters')
else
    vals = allVals(((sampleNum-1)*nParams+1):sampleNum*nParams);
end
fclose(fileID);
end