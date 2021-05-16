function [] = plot_channel_polygons(polygons,polyVertices2params, ...
                                    vertexCoords)
ftsize = 20;
color = 'k';
hAlignment = {'left','left','right','right'};
nhAlignment = numel(hAlignment);
vAlignment = {'top','bottom','bottom','top'};
nvAlignment = numel(vAlignment);
alignInd = 0;

for i = 1:numel(polygons)
    plot([polygons(i).XX(1,:),polygons(i).XX(1,1)],...
          [polygons(i).XX(2,:),polygons(i).XX(2,1)],'b-','linewidth',3)
    hold on
    text(mean(polygons(i).XX(1,:)),mean(polygons(i).XX(2,:)),num2str(i),...
         'fontsize',ftsize,'color','r')  
    
    for j = 1:size(polyVertices2params{i},2)      
        for k = 1:size(polyVertices2params{i},1)
            if (~isnan(polyVertices2params{i}(k,j)))
                alignInd = alignInd + 1;
                hAlign = hAlignment{rem(alignInd - 1,nhAlignment)+1};
                vAlign = vAlignment{rem(alignInd - 1,nvAlignment)+1};
                text(polygons(i).XX(1,j),polygons(i).XX(2,j), ...
                     ['[',num2str(polyVertices2params{i}(k,j)),']'],...
                    'color',color,'fontsize',ftsize,...
                    'VerticalAlignment', vAlign,...
                    'HorizontalAlignment',hAlign)
            end
        end
    end
    
end

color = [0,0.5,0];
for i = 1:size(vertexCoords,1)
    plot(vertexCoords(i,1),vertexCoords(i,2),...
         'color',color,'marker','o','markerfacecolor',color,...
         'markersize',8)
    text(vertexCoords(i,1),vertexCoords(i,2),num2str(i),...
        'color', color,'fontsize',ftsize)
end
xlabel('x', 'fontsize',ftsize);
ylabel('y', 'fontsize',ftsize);
axis image
set(gca,'fontsize',ftsize)
box on

end