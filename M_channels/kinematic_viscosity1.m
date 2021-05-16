%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Marcus Tan on 2/21/2014
%%% Copyright 2015 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% INPUT:
%%%     density: kg/m^3
%%%     T: temperature in K
%%%     model: model number
function [nu,DTnu] = kinematic_viscosity(density,T,model)
if model == 1
    factor = 0.0069/density;
    exponent = -8.3;
    Toffset = 273.15;
    nu = factor*(T/Toffset + 1.0)^exponent;
    if (nargout > 1)
        DTnu = (exponent*factor/Toffset)*(T/Toffset + 1.0)^(exponent-1);
    end
else
    error('unknown model')
end
end