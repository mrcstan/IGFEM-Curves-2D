function [] = plot_channel_polygons(polygons,vertexCoords,...
                                    nurbs,designParams, labelsOn)

ftsize = 40;
color = 'k';
hAlignment = {'left','left','right','right'};
nhAlignment = numel(hAlignment);
vAlignment = {'top','bottom','bottom','top'};
nvAlignment = numel(vAlignment);
alignInd = 0;

if (isempty(designParams))
    designParams.nParams = 0;
end

for i = 1:numel(polygons)
    XX = vertexCoords(polygons(i).connectivity,:);
    plot([XX(:,1);XX(1,1)],...
         [XX(:,2);XX(1,2)],'k-','linewidth',3)
    hold on
    %
    if (labelsOn)
        text(mean(XX(:,1)),mean(XX(:,2)),num2str(i),...
            'fontsize',ftsize,'color','r')
    end
    %
    %{
    designParamNum = designParams.vertices2params(polygons(i).connectivity,:);
    for j = 1:size(designParamNum,1)      
        for k = 1:size(designParamNum,2)
            if (~isnan(designParamNum(j,k)))
                alignInd = alignInd + 1;
                hAlign = hAlignment{rem(alignInd - 1,nhAlignment)+1};
                vAlign = vAlignment{rem(alignInd - 1,nvAlignment)+1};
                if (labelsOn)
                    text(XX(j,1),XX(j,2), ...
                         ['[',num2str(designParamNum(j,k)),']'],...
                        'color',color,'fontsize',ftsize,...
                        'VerticalAlignment', vAlign,...
                        'HorizontalAlignment',hAlign)
                end
            end
        end
    end
    %}
end

tol = 1e-10;
xmin = min(vertexCoords(:,1));
xmax = max(vertexCoords(:,1));
indLeft = abs((vertexCoords(:,1) - xmin)) < tol;
indRight = abs((vertexCoords(:,1) - xmax)) < tol;
ymin = min(vertexCoords(:,2));
ymax = max(vertexCoords(:,2));
indBottom = abs((vertexCoords(:,2) - ymin)) < tol;
indTop = abs((vertexCoords(:,2) - ymax)) < tol;

color = [0,0.5,0];
for i = 1:size(vertexCoords,1)
    %
    if ((indLeft(i) && ~indTop(i)) || (indRight(i) && ~indBottom(i)))
        plot(vertexCoords(i,1),vertexCoords(i,2),...
             'color',color,'marker','o','markerfacecolor',color,...
             'markersize',20)
    end
    %
    if ((indLeft(i) && indTop(i)) || (indRight(i) && indBottom(i)))
        plot(vertexCoords(i,1),vertexCoords(i,2),...
             'color','b','marker','d','markerfacecolor','b',...
             'markersize',20)
    end
    %
    if (labelsOn)
        text(vertexCoords(i,1),vertexCoords(i,2),num2str(i),...
            'color', 'k','fontsize',ftsize, ...
            'VerticalAlignment','bottom',...
            'HorizontalAlignment','left')
    end
    %
end

color = 'b';
marker = 'd';
%
for i = 1:designParams.nParams
    if (strcmpi(designParams.type{i},'CTRL_PT'))
        for j = 1:numel(designParams.channelNum{i})
            chanNum = designParams.channelNum{i}(j);
            cptNum = designParams.ctrlPtNum{i}(j);
            xx = nurbs(chanNum).coefs(1,cptNum)./nurbs(chanNum).coefs(4,cptNum);
            yy = nurbs(chanNum).coefs(2,cptNum)./nurbs(chanNum).coefs(4,cptNum);
            plot(xx,yy,'marker',marker,'markersize',20,'markerfacecolor',color);
            %{           
            alignInd = alignInd + 1;

            hAlign = hAlignment{rem(alignInd - 1,nhAlignment)+1};
            vAlign = vAlignment{rem(alignInd - 1,nvAlignment)+1};
            
            if (labelsOn)
                text(xx,yy,...
                    ['[',num2str(i),',',num2str(chanNum),',',num2str(cptNum),']'],...
                    'color',color,'fontsize',ftsize,...
                    'VerticalAlignment', vAlign,...
                    'HorizontalAlignment',hAlign)
            end
            %}
        end
        
    end
end
%
xlabel('x', 'fontsize',ftsize);
ylabel('y', 'fontsize',ftsize);
axis image
set(gca,'fontsize',ftsize)
axis off
delx = xmax - xmin;
dely = ymax - ymin;
xi = xmin - 0.1*delx;
yi = ymin - 0.1*dely;
delx = 1.2*delx;
dely = 1.2*dely;
rectangle('Position',[xi,yi,delx,dely],'edgeColor','k','clipping','off')
set(gca,'units','normalized','position',[0.15,0.15,0.7,0.7])
set(gcf,'color',[0.95,0.95,0.95])
end