% this function calculates volume fraction of microchannels assuming that
% they have circular cross sections
function vf =channel_volume(lengths,diameters,heights,crossSection)
if (strcmpi(crossSection,'circular'))
    vf = pi/4*sum(lengths.*diameters.^2);
elseif (strcmpi(crossSection,'square'))
    vf = sum(lengths.*diameters.^2);
elseif (strcmpi(crossSection,'rectangular'))
    vf = sum(lengths.*diameters.*heights);
else
    error('unknown cross section')
end
end