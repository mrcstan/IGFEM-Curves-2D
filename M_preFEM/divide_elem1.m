% Taken from the paper "Short Bisection Implementation in MATLAB" by Long
% Chen
% Modified by Marcus Tan on June 28, 2014
% INPUT:
%   markers: a length 3 array of the new node numbers on edge (2,3), (1,2)
%            and (1,3). If there is no new node an an edge, the marker is
%            set to 0
% OUTPUT:
%   modifiedElems: the number of modified as well as new elements
function [elem_nodes,modifiedElems] = divide_elem(elNum,elem_nodes,markers)
    modifiedElems = zeros(4,1);   
    p = [elem_nodes(elNum,:), markers(1)];
    nElem = size(elem_nodes,1);
    elem_nodes([elNum,nElem+1],:) = split(p);
    modifiedElems(1) = elNum;
    modifiedElems(2) = nElem+1;
    % divide right element further
    if (markers(2))
        p = [p(4),p(3),p(1),markers(2)];
        nElem = size(elem_nodes,1);
        elem_nodes([nElem,nElem+1],:) = split(p);
        modifiedElems(3) = nElem+1;
    end
    % divide left element further
    if (markers(3))
        p = [p(4),p(1),p(2),markers(3)];
        nElem = size(elem_nodes,1);
        elem_nodes([elNum,nElem+1],:) = split(p);
        modifiedElems(4) = nElem+1;
    end
    modifiedElems = modifiedElems(modifiedElems > 0);
end

function leftRightElem_nodes = split(p)
    leftRightElem_nodes(1,:) = [p(4),p(1),p(2)];
    leftRightElem_nodes(2,:) = [p(4),p(3),p(1)];
end