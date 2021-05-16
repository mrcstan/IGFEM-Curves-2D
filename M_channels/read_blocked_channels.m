%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Marcus Tan on 10/11/2015
%%% Copyright 2015 University of Illinois at Urbana-Champaign. 
%%% All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% INPUT: 
%%%     inputFile

%%% OUTPUT:
%%%     blockedSets: a cell array, of which each entry corresponds to a
%%%                  vector of blocked channel numbers
%%% FORMAT:
%{
<set 1 channel numbers>
<set 2 channel numbers>
...
<set n channel numbers>
%}
%%% FILE FORMAT EXAMPLE:
%{
2
3 
2 3 
%}
%%% 
function blockedSets = read_blocked_channels(inputFile)

if (~ischar(inputFile))
    error('input file name must be a character array');
end
fileID = fopen(inputFile,'r');
count = 0;
while ~feof(fileID)
    fgetl(fileID);
    count = count + 1;
end
fclose(fileID);

blockedSets = cell(count,1);
count = 0;
fileID = fopen(inputFile,'r');
while ~feof(fileID)
    line = fgetl(fileID);
    count = count + 1;
    blockedSets{count} = str2num(line);
end
 
fclose(fileID); 
end