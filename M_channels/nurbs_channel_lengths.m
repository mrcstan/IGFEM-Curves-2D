%%% Created by Marcus Tan in 2014
%%% Modified on 1/7/2015
%%% Copyright 2014 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function calculates the length of channels described by NURBS
% INPUT: 
%   channels: structure array of nurbs
%   ngpt: number of gauss points
% OUTPUT: (last 3 outputs are )
%   lengths:
%   lineSegs: for intersect_edges function, only valid for linear NURBS
%   segLengths: length of each line segment, only valid for linear NURBS
function [lengths,lineSegs,segLengths]...
            = nurbs_channel_lengths(channels,ngpts)
if (~isstruct(channels))
    error('input must be a nurbs structure array')
end
nChannels = numel(channels.nurbs);
lengths = zeros(nChannels,1);
if (nargout > 1)
    lineSegs = cell(nChannels,1);
    segLengths = cell(nChannels,1);
    for i = 1: nChannels
        if (channels.nurbs(i).order == 2) 
            segLengths{i} = sqrt(sum((channels.nurbs(i).coefs(1:2,2:end)...
                             -channels.nurbs(i).coefs(1:2,1:end-1)).^2,1));
            lengths(i) = sum(segLengths{i});
            lineSegs{i} = [channels.nurbs(i).coefs(1:2,1:end-1);...
                           channels.nurbs(i).coefs(1:2,2:end)]';  

        else
            lengths(i) = nurbs_arc_length(ngpts,channels.nurbs(i));
        end   
    end
else
   for i = 1: nChannels
        if (channels.nurbs(i).order == 2)        
            lengths(i) = sum(sqrt(sum((channels.nurbs(i).coefs(1:2,2:end)...
                             -channels.nurbs(i).coefs(1:2,1:end-1)).^2,1)));

        else
            lengths(i) = nurbs_arc_length(ngpts,channels.nurbs(i));
        end   
   end
end
end