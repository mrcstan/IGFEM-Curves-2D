%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Marcus Tan on 9/17/2013
%%% Last modified date: 9/17/2013
%%% Copyright 2013 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this function construct convex hulls of nurbs curves. when points are collinear
% a rectangular container of half width=halfWidth is constructed instead
% INPUT:
%   flag:  1: the convexhull consists of a series of convex hulls
%          2: the convexhull consists of the union of the smaller convex hulls 
function itrface=construct_nurbs_convexhull(itrface,hullParams)

for j=1:itrface.nNurbs
    ctrlPt=[itrface.nurbs(j).coefs(1,:);
            itrface.nurbs(j).coefs(2,:)];
    % if all control points are collinear, just construct a straight line    
    if(iscollinear(ctrlPt,hullParams.collinearTol))
        itrface.nurbs(j).nConvHull=1;
        itrface.nurbs(j).iscollinear{1} = true;
        itrface.nurbs(j).hull{1} = line_segment_of_col_pts(ctrlPt);  
    else
        itrface.nurbs(j).nConvHull=itrface.nurbs(j).number-itrface.nurbs(j).order+1;
        for k = 1:itrface.nurbs(j).nConvHull
            m=k+itrface.nurbs(j).order-1;
            ctrlPt=[itrface.nurbs(j).coefs(1,k:m)./itrface.nurbs(j).coefs(4,k:m);
                    itrface.nurbs(j).coefs(2,k:m)./itrface.nurbs(j).coefs(4,k:m)];
            
            itrface.nurbs(j).iscollinear{k} = iscollinear(ctrlPt,hullParams.collinearTol);    
            % if control points are collinear, construct a rectangular hull         
            if(itrface.nurbs(j).iscollinear{k})          
                
                %itrface.nurbs(j).hull{k} = line_segment_of_col_pts(ctrlPt);
                %itrface.nurbs(j).hull{k} = reshape(itrface.nurbs(j).hull{k},2,2);
                itrface.nurbs(j).hull{k} = construct_rect_for_col_pt(ctrlPt,hullParams.halfWidth);
            else
                %itrface.nurbs(j).hull{k} =ctrlPt;
                itrface.nurbs(j).hull{k} = scale_polygon(ctrlPt,hullParams.expandFactor);
            end
        end
        
        if (hullParams.type == 2)
            X = itrface.nurbs(j).hull{1}(1,:);
            Y = itrface.nurbs(j).hull{1}(2,:);
            [X,Y] = poly2cw(X,Y);
            for k = 2:itrface.nurbs(j).nConvHull
                %{
                if (~itrface.nurbs(j).iscollinear{k})
                    
                else
                    X2 = [itrface.nurbs(j).hull{k}(1),itrface.nurbs(j).hull{k}(3)];
                    Y2 = [itrface.nurbs(j).hull{k}(2),itrface.nurbs(j).hull{k}(4)];
                end
                %}
                [X2,Y2] = poly2cw(itrface.nurbs(j).hull{k}(1,:),itrface.nurbs(j).hull{k}(2,:));
                [X,Y] = polybool('union',X,Y,X2,Y2);
            end
            itrface.nurbs(j).hull{1} = [X;Y];
            itrface.nurbs(j).nConvHull = 1;
            itrface.nurbs(j).iscollinear{1} =  false;
        end                            
    end
end
end

function line = line_segment_of_col_pts(pt)
npt=size(pt,2);
if(npt < 2)
    error('construct_rect_for_col_pt: must have at least two points');
end
dirVec = pt(1:2,2)-pt(1:2,1);
normDirVec = norm(dirVec,2);
if(normDirVec == 0)
    error('end_pts_of_col_pts: some points are coincident')
else
    dirVec = dirVec/normDirVec;
end    

%sort the points along the collinear line such that in ascending order
%of distance from an end point
signedDist = zeros(1,npt);    
for i=2:npt
    signedDist(i)=dot(pt(1:2,i)-pt(1:2,1),dirVec);
end
[~,ind]=sort(signedDist);
line = [pt(1:2,ind(1))',pt(1:2,ind(end))'];

end
