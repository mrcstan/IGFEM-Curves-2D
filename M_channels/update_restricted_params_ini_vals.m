%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Marcus Tan on 1/19/2015
%%% Copyright 2015 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this function equates the restricted parameters initial values to their
% corresponding design parameters. 
% in restrictedParams.paramPairs, first column is the designParam
% numbers and second column is the restricted param numbers
function restrictedParamIniVals ...
    = update_restricted_params_ini_vals(designParams, ...
                                        restrictedParams)

restrictedParamIniVals = restrictedParams.iniVals;
for i = 1:size(restrictedParams.paramPairs,1)
    restrictedParamIniVals(restrictedParams.paramPairs(i,2)-designParams.nParams) ...
        = designParams.iniVals(restrictedParams.paramPairs(i,1)); 
end                                    
end                                                           