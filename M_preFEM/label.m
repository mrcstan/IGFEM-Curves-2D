% Taken from the paper "Short Bisection Implementation in MATLAB" by Long
% Chen
% Modified by Marcus Tan on June 22, 2014
function elem=label(node,elem)

    edgelength(:,1) = (node(elem(:,3),1)-node(elem(:,2),1)).^2 ...
                    + (node(elem(:,3),2)-node(elem(:,2),2)).^2;
    edgelength(:,2) = (node(elem(:,1),1)-node(elem(:,3),1)).^2 ...
                    + (node(elem(:,1),2)-node(elem(:,3),2)).^2;
    edgelength(:,3) = (node(elem(:,2),1)-node(elem(:,1),1)).^2 ...
                    + (node(elem(:,2),2)-node(elem(:,1),2)).^2;            
          
    [~,I] = max(edgelength,[],2);
    elem((I==2),[1 2 3]) = elem((I==2),[2 3 1]);
    elem((I==3),[1 2 3]) = elem((I==3),[3 1 2]);
end
