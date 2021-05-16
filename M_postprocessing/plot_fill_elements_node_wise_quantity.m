%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Marcus Tan on 6/12/2013
%%% Last modified date: 6/12/2013
%%% Copyright 2013 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_fill_elements_node_wise_quantity(nodeCoords,elem,UUR,scale)
nscale=4; % the position of the length scale in the arguments
if(nargin<nscale)
    scale.length=1;
    scale.uo=1;
    scale.ud=1;
end

figure

ftsize=30;
%subd=30;
specs.curve = 'w-';
specs.width = 0.5;
hold on;

for i = 1: elem.n_elem
    if (elem.parent(i).type > 1)
        for j=1:numel(elem.parent(i).child)
            X=nodeCoords(elem.parent(i).child(j).nodes,1)/scale.length;
            Y=nodeCoords(elem.parent(i).child(j).nodes,2)/scale.length;
            hullInd=convhull(X,Y);
            p = patch(X(hullInd),Y(hullInd),(UUR(elem.parent(i).child(j).nodes(hullInd))-scale.uo)/scale.ud);
            set(p,'EdgeColor','none')
        end
        
        %{
        for j=1:numel(elem.parent(i).channelNum)  
            elem.parent(i).nurbsSeg(j).coefs(1:2,:) = ...
                elem.parent(i).nurbsSeg(j).coefs(1:2,:)/scale.length;
            nurbs_plot(elem.parent(i).nurbsSeg(j),subd,specs);
        end
        %}
    else
         X=nodeCoords(elem.elem_node(i,:),1)/scale.length; 
         Y=nodeCoords(elem.elem_node(i,:),2)/scale.length;     
         p = patch(X, Y, (UUR(elem.elem_node(i,:))-scale.uo)/scale.ud); 
         set(p,'EdgeColor','none')
         if(elem.parent(i).type == 1)
             for j=1:numel(elem.parent(i).channelNum)
                 %{
                 elem.parent(i).nurbsSeg(j).coefs(1:2,:) = ...
                    elem.parent(i).nurbsSeg(j).coefs(1:2,:)/scale.length;
                nurbs_plot(elem.parent(i).nurbsSeg(j),subd,specs);
                 %}
                 X = nodeCoords(elem.parent(i).lineSource(1:2,j),1)/scale.length;
                 Y = nodeCoords(elem.parent(i).lineSource(1:2,j),2)/scale.length;
                 plot(X,Y,specs.curve,'linewidth',specs.width')
             end
         end
    end
    
end
%title('Temperature distribution','fontsize',ftsize)
xlabel('x','fontsize',ftsize);
ylabel('y','fontsize',ftsize,'Rotation',0);
axis xy;
axis image;
axis tight;
%caxis([0,0.35])
h = colorbar;
ylabel(h,'T^*','fontsize',ftsize,'Rotation',0)

set(gca,'fontsize',ftsize)

end