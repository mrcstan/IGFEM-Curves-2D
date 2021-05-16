%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Marcus Tan on 2/21/2014
%%% Copyright 2015 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% INPUT:
%%%     nu1: kinematic viscosity used only for default model
%%%     density: kg/m^3
%%%     T: temperature in K
%%%     model: 
%%%         1: correlation of the kinematic viscosity with temperature
%%%         0: constant kinematic viscosity given by nu1
function [nu,DTnu] = kinematic_viscosity(nu1,density,T,model)
switch model
    case 0
        nu = nu1;
        DTnu = 0.0;
    case 1
        factor = 0.0069/density;
        exponent = -8.3;
        Toffset = 273.15;
        nu = factor*(T/Toffset + 1.0)^exponent;
        if (nargout > 1)
            DTnu = (exponent*factor/Toffset)*(T/Toffset + 1.0)^(exponent-1);
        end
    otherwise
        error('unknown viscosity model')
end
end