%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Marcus Tan on 11/6/2014
%%% Last modified date: 11/6/2014
%%% Copyright 2014 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this function generates size(bounds,1) random numbers between bounds(:,1)
% and bounds(:,1) inclusive
function randVals = randomized_bounded_values(bounds)
rng('shuffle')
randVals = bounds(:,1) + (bounds(:,2)-bounds(:,1)).*rand(size(bounds,1),1);
end