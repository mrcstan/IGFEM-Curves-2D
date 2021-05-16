%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Marcus Tan on 8/5/2013
%%% Last modified date: 12/27/2013
%%% Copyright 2013 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% INPUT: 
function plot_fill_elements_elem_wise_error_or_norm(nodeCoords,elem,normType,plotType)

if (strcmpi(normType,'L2'))
    normID = 1;
elseif (strcmpi(normType,'H1'))
    normID = 2;
else
    error('unknown norm type')
end

if (strcmpi(plotType,'error'))
    plotID = 1;
elseif (strcmpi(plotType,'norm'))
    plotID = 2;
else
    error('unknown plot type')
end
edge_color = 'none';
ftsize = 24;
specs.curve = 'k-';
specs.width = 2;
subd = 30;


figure
hold on;
for i = 1: elem.n_elem        
    % NURBS IGFEM
     if (elem.parent(i).type == 3)  
         for c=1:numel(elem.parent(i).child)
            X=nodeCoords(elem.parent(i).child(c).nodes,1);
            Y=nodeCoords(elem.parent(i).child(c).nodes,2);
            hullInd=convhull(X,Y);
            if (plotID == 1)
                if (normID == 1)
                    p = patch(X(hullInd),Y(hullInd),elem.parent(i).child(c).L2);
                else
                    p = patch(X(hullInd),Y(hullInd),elem.parent(i).child(c).H1);
                end
            else
                if (normID == 1)
                    p = patch(X(hullInd),Y(hullInd),elem.parent(i).child(c).normL2);
                else
                    p = patch(X(hullInd),Y(hullInd),elem.parent(i).child(c).normH1);
                end
            end
            set(p,'EdgeColor',edge_color);
         end
         
         for j=1:numel(elem.parent(i).channelNum)  
            nurbs_plot(elem.parent(i).nurbsSeg(j),subd,specs);
         end
     elseif (elem.parent(i).type == 2)
          X=nodeCoords(elem.parent(i).nodes(1:3),1);
         Y=nodeCoords(elem.parent(i).nodes(1:3),2);
         if (plotID == 1)
             if (normID == 1)
                p = patch(X,Y,elem.parent(i).L2);
             else
                p = patch(X,Y,elem.parent(i).H1);
             end
         else
             if (normID == 1)
                p = patch(X,Y,elem.parent(i).normL2);
             else
                p = patch(X,Y,elem.parent(i).normH1);
             end
         end
         set(p,'EdgeColor',edge_color);
     else
         X = nodeCoords(elem.elem_node(i,:),1); 
         Y = nodeCoords(elem.elem_node(i,:),2);
         if(plotID==1)
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
if(plotID==1)
    if (normID == 1)
        title('L^2 error');
    else
        title('H^1 error');
    end
else
    if (normID == 1)
        title('L^2 norm of analytical solution');
    else
        title('H^1 norm of analytical solution');
    end
end
axis xy;
axis image;
axis tight;
colorbar;
set(gca,'fontsize',ftsize)


end