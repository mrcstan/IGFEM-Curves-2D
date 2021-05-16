%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Marcus Tan on 1/16/2015
%%% Copyright 2014 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% merge 2 triangles to form a quadrilatel if 
% (i) the triangles do not share a constraint edge.
% (ii) the resulting quadrilateral do not have an interior angle> 180 deg
% INPUT:
%   TR: connectivitiy of triangles, each row corresponding to 1 triangle
%   pts: coordinates of the vertices, each row corresponding to a vertex
%   constraints: constrained edges, each col correponding to an edge
% OUTPUT:
%   newTR: connectivity of unmerged triangles
%   rec: connectivity of quadrilaterals, each row corresponding to 1
%        quadrilateral
function [newTR,rec] = merge_triangles(TR,pts,constraints)
nTri = size(TR,1);
delTR = false(nTri,1);
rec = nan(nTri,4);
nRec = 0;

tol = 1e-8;
for i=1:(nTri-1)
    for j=i+1:nTri
        % check if the triangle has already been merged.
        % if so skip checking this triangle
        if(~delTR(i) && ~delTR(j))
        %if(isempty(find(delTR==j, 1)))
            cij=intersect(TR(i,:),TR(j,:));
            % if two triangles share two common nodes, 
            % then they must share an edge
            if(length(cij)==2)
                edgeConstraint=false;
                % check if the edge is a constraint edge    
                for k = 1:size(constraints,2)
                    com=setdiff(cij, constraints(:,k));
                    if(isempty(com))
                        edgeConstraint=true;
                        break;
                    end
                end
     
                if(~edgeConstraint)
                    [dij, di]=setdiff(TR(i,:),TR(j,:));
                    [dji, dj]=setdiff(TR(j,:),TR(i,:));
                    % check if resulting rectangle will have interior 
                    % angle greater than 180 deg
                    vec1 = pts(dij,:)-pts(cij(1),:);
                    vec2 = pts(cij(2),:)-pts(cij(1),:);
                    vec3 = pts(dji,:)-pts(cij(1),:);
                    r1=vec1*vec2'/(norm(vec1,2)*norm(vec2,2));
                    r2=vec3*vec2'/(norm(vec3,2)*norm(vec2,2));
                    sinAngle=sqrt(1-r1^2)*r2+sqrt(1-r2^2)*r1;
                    if(sinAngle<=tol)
                        break
                    end
                    vec1 = pts(dij,:)-pts(cij(2),:);
                    vec2 = pts(cij(1),:)-pts(cij(2),:);
                    vec3 = pts(dji,:)-pts(cij(2),:);
                    r1 = vec1*vec2'/(norm(vec1,2)*norm(vec2,2));
                    r2 = vec3*vec2'/(norm(vec3,2)*norm(vec2,2));
                    sinAngle = sqrt(1-r1^2)*r2+sqrt(1-r2^2)*r1;
                    if(sinAngle<=tol)
                        break
                    end
                    dir=di+1;
                    if(dir>3)
                        dir=1;
                    end
                    djr=dj+1;
                    if(djr>3)
                        djr=1;
                    end
                    nRec=nRec+1;
                    rec(nRec,:)=[dij(1),TR(i,dir),dji(1),TR(j,djr)]; 
                    %delTR=[delTR,i,j];
                    delTR(i) = true;
                    delTR(j) = true;
                    break; 
                end 
            end
        end
    end
end
newTR = TR(~delTR,:);

rec((nRec+1):end,:) = [];
end