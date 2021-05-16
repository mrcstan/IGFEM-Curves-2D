%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Marcus Tan on 1/22/2014
%%% Last modified date: 10/22/2014
%%% Copyright 2014 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check whether a real number is in the sorted real number list. if yes,
% return the location of the number
% if not, insert the real number in the appropriate location to maintain the
% sorted list. The location of the smallest number in the original list
% greater than equal the query number is returned.
function [inList,loc,list] = find_and_update_sorted_real(query,list,tol)
if(isempty(list))
    inList = false;
    loc = 1;
    list = query;
else
    nParams = numel(list);
    %{
    if(abs(list(1)-query)<tol)
        inList = true;
        loc = 1;        
        return
    elseif(query < list(1))
        inList = false;
        loc = 1;        
        list = [query,list];
        return
    end
    for i = 2:nParams
        if(abs(list(i)-query)<tol)
            inList = true;
            loc = i;
            return
        elseif(query < list(i))
            inList = false;
            loc = i;
            list = [list(1:i-1),query,list(i:end)];
            return
        end
    end
    inList = false;
    loc = nParams+1;
    list(nParams+1) = query;
    %}
    [loc,inList] = binary_search(query,list,1,nParams);
    if ~inList
        if abs(query-list(loc)) < tol
            % query is slightly smaller than list(loc)
            inList = true;
            return
        end
        if loc > 1  && abs(query-list(loc-1)) < tol
            % query is slightly larger than list(loc-1)
            inList = true;
            loc = loc-1;
            return
        end    
        if loc == nParams && query > list(loc)
            % handle the case where the query is greater than any entry of the list
            loc = nParams+1; 
            list(nParams+1) = query;
        else
            list = [list(1:loc-1),query,list(loc:end)];
        end
    end
    
end