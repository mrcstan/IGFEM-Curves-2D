%%% Created by Marcus Tan on 9/21/2014
%%% Last modified date: 9/21/2014
%%% Copyright 2014 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% binary search
function [imin,found] = binary_search(key,sortedA,imin,imax)

    while (imin < imax)
        imid = imin + floor((imax-imin)/2); % prevents overflow
        
        if (sortedA(imid) < key)
            imin = imid + 1;
        else
            imax = imid;
        end
    end
    
    
    if (sortedA(imin) == key)
        found = true;
    else
        found = false;
    end
end
