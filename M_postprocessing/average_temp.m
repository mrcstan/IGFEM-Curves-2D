 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Marcus Tan on 2/18/2014
%%% Modified by Marcus Tan
%%% Last modified date: 2/18/2014
%%% Copyright 2014 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculates an approximate area weighted average temperature
function ave = average_temp(elem,nodeCoords,UUR)

nanInd = isnan(UUR);
if (any(nanInd))
	warning('removing NAN values from UUR') 
	UUR(nanInd) = 0.0; % set nan to 0 in case the equation is singular
end

totArea = 0.0;
ave = 0.0;

for i = 1:elem.n_elem
    if(elem.parent(i).type >1)
        for c = 1:numel(elem.parent(i).child)
            childArea = abs(polygonArea(nodeCoords(elem.parent(i).child(c).nodes,:)));
            totArea = totArea + childArea;
            ave = ave+mean(UUR(elem.parent(i).child(c).nodes))*childArea;
        end
    else
        elemArea = abs(polygonArea(nodeCoords(elem.elem_node(i,:),:)));
        totArea = totArea + elemArea;
        ave = ave+mean(UUR(elem.elem_node(i,:)))*elemArea;
    end
end
ave = ave/totArea;
end