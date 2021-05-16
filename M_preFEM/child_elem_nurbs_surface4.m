%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Marcus Tan on 8/25/13
%%% Modified by Marcus Tan
%%% Last modified date: 7/15/14
%%% Copyright 2013 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this function calculates the nurbs surface representing a triangular
% element or quadrilateral element with one or two edges described by a
% nurbs curve. 
% ASSUMPTION: i) The case where more than two edges are nurbs 
%               curves are not considered.
%             ii) When there is branching, the maximum control points per 
%                 NURBS segment is 3. This is because only 1 interior point
%                 is created by this code. 
% INPUT:
%   xel: a 2x3(triangular) or 2x4(quadrilateral) array of the coordinates
%   of the vertices of the element.
%   lineSourceLoc: a length 2 (1 nurbs curve) or 4 (2 nurbs curves) of the 
%               local vertex numbers of end points of the nurbs segment
%              Ex. [1,2,4,3] means there are two nurbs curves. The first
%              curve points from 1 to 2, the second curve points from 4 to
%              3.
%   nurbsSeg: a single nurbs data structure or a length 2 array of nurbs
%   data structure. Two nurbs curves must be in the same "direction".
% OUTPUT:
%   nurbsSurf
function [nurbsSurf,uv,uvParam]=child_elem_nurbs_surface(xel,lineSourceLoc,nurbsSeg)
if(size(xel,2)==3)
    shape=1;
elseif(size(xel,2)==4)
    shape=2;
end
nNurbsSeg = length(nurbsSeg);
if(nNurbsSeg>2)
    error('child_elem_nurbs_surface: case of child element with more than 3 line sources not treated')
end
if(shape==1)
    if(nNurbsSeg == 1)
        otherLocNode=setdiff(1:3,lineSourceLoc(1,1:2));
        if(isCCW(lineSourceLoc(1,1:2),3))
            surfPt(:,:,1)=nurbsSeg(1).coefs;
        else
            surfPt(:,:,1)=nurbsSeg(1).coefs(:,end:-1:1);
        end
        surfPt(:,:,2)=repmat([xel(1,otherLocNode);xel(2,otherLocNode);0;1],...
                             [1,nurbsSeg(1).number]);
        uknot=nurbsSeg(1).knots;
        vknot=[0,0,1,1];
        uv = 2;
        uvParam = 0;
    elseif(nNurbsSeg == 2)
        % check if the two NURBS segment have the same number of control
        % points
        if (nurbsSeg(1).number ~= nurbsSeg(2).number)
            error('2 NURBS segments do not have the same number of control points')
        end
        if (nurbsSeg(1).order ~= nurbsSeg(2).order)
            error('2 NURBS segments do not have the same order')
        end
        if(isCCW(lineSourceLoc(1,1:2),3))
            surfPt(:,:,1)=nurbsSeg(1).coefs;
        else
            surfPt(:,:,1)=nurbsSeg(1).coefs(:,end:-1:1);
        end
        if(isCCW(lineSourceLoc(2,1:2),3))    
            surfPt(:,:,2)=nurbsSeg(2).coefs(:,end:-1:1); 
        else         
            surfPt(:,:,2)=nurbsSeg(2).coefs;
        end       
        uknot=nurbsSeg(1).knots;
        vknot=[0,0,1,1];
        uv = [2;2];
        uvParam = [0;1];
    elseif(nNurbsSeg == 0)
        surfPt(1:2,1,1) = xel(1:2,1);
        surfPt(1:2,2,1) = xel(1:2,2);
        surfPt(1:2,1,2) = xel(1:2,3);
        surfPt(1:2,2,2) = xel(1:2,3);
        uknot = [0,0,1,1];
        vknot = [0,0,1,1];
        uv = [];
        uvParam = [];
    end
elseif(shape == 2)
    if(nNurbsSeg == 1)
        if(isCCW(lineSourceLoc(1,1:2),4))
            surfPt(:,:,1) = nurbsSeg(1).coefs;
        else
            surfPt(:,:,1) = nurbsSeg(1).coefs(:,end:-1:1);
        end
        nonSourceLocNode = setdiff(1:4,lineSourceLoc(1,1:2));
        if(isCCW(nonSourceLocNode,4))
            nonSourceLocNode = nonSourceLocNode(2:-1:1);
        end
        x = linspace(xel(1,nonSourceLocNode(1)),xel(1,nonSourceLocNode(2)),nurbsSeg(1).number);
        y = linspace(xel(2,nonSourceLocNode(1)),xel(2,nonSourceLocNode(2)),nurbsSeg(1).number);
        z = zeros(1,nurbsSeg(1).number);
        w = ones(1,nurbsSeg(1).number);
        surfPt(:,:,2) = [x;y;z;w];
        uknot = nurbsSeg(1).knots;
        vknot = [0,0,1,1];
        uv = 2;
        uvParam = 0;
    elseif(nNurbsSeg == 2)
        error('cannot handle quadrilateral child elements with 2 line sources') 
        %{
        % the two nurbs segmets do not meet each other.
        commonLoc = intersect(lineSourceLoc(1,1:2),lineSourceLoc(3:4),'stable');
        if(isempty(commonLoc))
            if(isCCW(lineSourceLoc(1,1:2),4))
                surfPt(:,:,1) = nurbsSeg(1).coefs;
            else
                surfPt(:,:,1) = nurbsSeg(1).coefs(:,end:-1:1);
            end
            if(isCCW(lineSourceLoc(2,1:2),4))
                surfPt(:,:,2) = nurbsSeg(2).coefs(:,end:-1:1);            
            else
                surfPt(:,:,2) = nurbsSeg(2).coefs;
            end        
            uknot = nurbsSeg(1).knots;
            vknot = [0,0,1,1];
        % the two nurbs segments intersect at one of their end points    
        else
            otherNode = setdiff(1:4,lineSourceLoc(:));
            isCCW1 = isCCW(lineSourceLoc(1,1:2),4);
            isCCW2 = isCCW(lineSourceLoc(2,1:2),4);
            if(~isCCW1)
                lineSourceLoc(1,1:2)=lineSourceLoc(2:-1:1);
            end
            if(~isCCW2)     
                lineSourceLoc(2,1:2)=lineSourceLoc(2:-1:1);  
            end        
            
            commonLocOnSeg1 = find(lineSourceLoc(1,1:2)==commonLoc,1,'first'); 
            if(commonLocOnSeg1 == 1) 
                % swap line sources
                [lineSourceLoc(2,1:2),lineSourceLoc(1,1:2)] ...
                        = deal(lineSourceLoc(1,1:2),lineSourceLoc(2,1:2)); 
                [isCCW2,isCCW1] = deal(isCCW1,isCCW2);
                segNum1=2;
                segNum2=1;
                
            else
                segNum1=1;
                segNum2=2;             
            end
            sLocNode=lineSourceLoc(1);
            eLocNode=lineSourceLoc(4);
            
            surfPt=zeros(4,nurbsSeg(segNum1).number,nurbsSeg(segNum2).number);
            x=linspace(xel(1,otherNode),xel(1,eLocNode),nurbsSeg(segNum1).number);
            y=linspace(xel(2,otherNode),xel(2,eLocNode),nurbsSeg(segNum1).number);
            z=zeros(1,nurbsSeg(segNum1).number);
            w=ones(1,nurbsSeg(segNum1).number);
            surfPt(:,:,nurbsSeg(segNum2).number)=[x;y;z;w];
            if(isCCW1)
                surfPt(:,:,1)=nurbsSeg(segNum1).coefs;
            else
                surfPt(:,:,1)=nurbsSeg(segNum1).coefs(:,end:-1:1); 
            end
            uknot=nurbsSeg(segNum1).knots;
            
            x=linspace(xel(1,sLocNode),xel(1,otherNode),nurbsSeg(segNum2).number);
            y=linspace(xel(2,sLocNode),xel(2,otherNode),nurbsSeg(segNum2).number);
            z=zeros(1,nurbsSeg(segNum2).number);
            w=ones(1,nurbsSeg(segNum2).number);
            surfPt(:,1,:)=[x;y;z;w];
            if(isCCW2)
                surfPt(:,nurbsSeg(segNum1).number,:)=nurbsSeg(segNum2).coefs;
            else
                surfPt(:,nurbsSeg(segNum1).number,:)=nurbsSeg(segNum2).coefs(:,end:-1:1);
            end
            vknot=nurbsSeg(segNum2).knots;
            
            % define interior control points which are distributed as
            % evenly as possible. Currently only works for 3 control points
            % in both u and v directions.
            surfPt = create_interior_control_points(surfPt);
        end
        %}
    elseif(nNurbsSeg == 0)
        surfPt(1:2,1,1) = xel(1:2,1);
        surfPt(1:2,2,1) = xel(1:2,2);
        surfPt(1:2,2,2) = xel(1:2,3);
        surfPt(1:2,1,2) = xel(1:2,4);
        uknot = [0,0,1,1];
        vknot = [0,0,1,1];
        uv = [];
        uvParam = [];
    end
end
nurbsSurf=nrbmak(surfPt,{uknot,vknot});
end



%{
function surfPt = create_interior_control_points(surfPt)
    sz=size(surfPt);
    if(sz(2)==3 && sz(3)==3)
        surfPt(1,2,2)=sum(surfPt(1,:));
        surfPt(1,2,2)=surfPt(1,2,2)/8;
        surfPt(2,2,2)=sum(surfPt(2,:));
        surfPt(2,2,2)=surfPt(2,2,2)/8;
        surfPt(4,2,2)=1;
    end
end
%}