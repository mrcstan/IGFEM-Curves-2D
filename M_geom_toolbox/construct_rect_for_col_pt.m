%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Marcus Tan on 8/20/2013
%%% Last modified date: 1/8/2014
%%% Copyright 2013 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this function construct a rectangle of half width halfWidth containing
% all the collinear points
% INPUT: ain 
%   pt: 2xn array of collinear point coordinates
%   halfWidth
% OUTPUT:
%   rect: 2x4 array of the rectangle arranged in CCW or CW order
function rect=construct_rect_for_col_pt(pt,halfWidth)
    npt=size(pt,2);
    if(npt < 2)
        error('construct_rect_for_col_pt: must have at least two points');
    end
    dirVec=[pt(1,2)-pt(1,1);pt(2,2)-pt(2,1)];
    normDirVec = norm(dirVec,2);
    if(normDirVec == 0)
        error('construct_rect_for_col_pt: some points are coincident')
    else
        dirVec = dirVec/normDirVec;
    end    
    
    %sort the points along the collinear line such that in ascending order
    %of distance from an end point
    signedDist = zeros(1,npt);    
    for i=2:npt
        signedDist(i)=dot([pt(1,i)-pt(1,1);pt(2,i)-pt(2,1)],dirVec);
    end
    [~,ind]=sort(signedDist);
    pt(1,:)=pt(1,ind);
    pt(2,:)=pt(2,ind);
    % construct the rectangular container
    dirVec = [pt(1,end)-pt(1,1);pt(2,end)-pt(2,1)];
    dirVec=dirVec*halfWidth/norm(dirVec,2);
    normVec=[-dirVec(2);dirVec(1)];
    rect=zeros(2,4);
    rect(:,1)=[pt(1,1)-normVec(1)-dirVec(1);pt(2,1)-normVec(2)-dirVec(2)];
    rect(:,2)=[pt(1,1)+normVec(1)-dirVec(1);pt(2,1)+normVec(2)-dirVec(2)];
    rect(:,3)=[pt(1,end)+normVec(1)+dirVec(1);pt(2,end)+normVec(2)+dirVec(2)];
    rect(:,4)=[pt(1,end)-normVec(1)+dirVec(1);pt(2,end)-normVec(2)+dirVec(2)];
end
