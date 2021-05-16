%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Marcus Tan on 8/5/2013
%%% Last modified date: 9/18/2013
%%% Copyright 2013 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% INPUT: 
%   type: 1, element-wise normalized L2 and H1 error
%         2, element-wise L2 and H1 norms of solution
function plot_fill_elements_elem_wise_error_or_norm(allNodeCoords,elem,type,scale)


edge_color = 'none';
ftsize = 24;
specs.curve = 'k-';
specs.width = 2;
subd = 30;
if(nargin<4)
    scale.length=1.0;
end

figure
hold on;
for i = 1: elem.n_elem         
     if (elem.parent(i).type > 1 && type==1)              
         for c=1:numel(elem.parent(i).child)
            X=allNodeCoords(elem.parent(i).child(c).nodes,1)/scale.length;
            Y=allNodeCoords(elem.parent(i).child(c).nodes,2)/scale.length;
            hullInd=convhull(X,Y);
            p = patch(X(hullInd),Y(hullInd),elem.parent(i).child(c).L2);
            set(p,'EdgeColor',edge_color);
         end
         for j=1:numel(elem.parent(i).channelNum)  
            elem.parent(i).nurbsSeg(j).coefs(1:2,:) = ...
                elem.parent(i).nurbsSeg(j).coefs(1:2,:)/scale.length;
            nurbs_plot(elem.parent(i).nurbsSeg(j),subd,specs);
        end
     else
         X=allNodeCoords(elem.elem_node(i,:),1)/scale.length; 
         Y=allNodeCoords(elem.elem_node(i,:),2)/scale.length;
         if(type==1)
            p = patch(X,Y,elem.parent(i).L2);
            set(p,'EdgeColor',edge_color);
         else
            p = patch(X,Y,elem.parent(i).normL2,'EdgeColor', edge_color);
            set(p,'EdgeColor',edge_color);
         end
     end
end

xlabel('x','fontsize',ftsize);
ylabel('y','fontsize',ftsize);
if(type==1)
    title('L^2 error');
else
    title('L^2 norm of analytical solution');
end
axis xy;
axis image;
axis tight;
colorbar;
set(gca,'fontsize',ftsize)

figure 
hold on

for i = 1: elem.n_elem    
      if(elem.parent(i).type > 1 && type==1)           
        for c=1:numel(elem.parent(i).child)
            X=allNodeCoords(elem.parent(i).child(c).nodes,1)/scale.length;
            Y=allNodeCoords(elem.parent(i).child(c).nodes,2)/scale.length;
            hullInd=convhull(X,Y);
            p=patch(X(hullInd),Y(hullInd),elem.parent(i).H1);    
             set(p,'EdgeColor',edge_color);
        end
        for j=1:numel(elem.parent(i).channelNum)  
            %elem.parent(i).nurbsSeg(j).coefs(1:2,:) = ...
            %    elem.parent(i).nurbsSeg(j).coefs(1:2,:)/scale.length;
            nurbs_plot(elem.parent(i).nurbsSeg(j),subd,specs);
        end
     else
         X=allNodeCoords(elem.elem_node(i,:),1)/scale.length; 
         Y=allNodeCoords(elem.elem_node(i,:),2)/scale.length;
         if(type==1)
            p=patch(X,Y,elem.parent(i).H1);
           set(p,'EdgeColor',edge_color);
         else
            p=patch(X,Y,elem.parent(i).normH1);
            set(p,'EdgeColor',edge_color);
         end
     end
end




xlabel('x','fontsize',ftsize);
ylabel('y','fontsize',ftsize);
if(type==1)
    title('H^1 error');
else
    title('H^1 norm of analytical solution');
end
axis xy;
axis image;
axis tight;
colorbar;
set(gca,'fontsize',ftsize)

end