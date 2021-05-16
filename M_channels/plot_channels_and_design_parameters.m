%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Marcus Tan on 7/10/2014
%%% Last modified date: 7/10/2014
%%% Copyright 2014 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this function plots the channel network if the nodal positions 
% or nurbs descriptions of the channels are given
% INPUT:
%   chan_nodes: a number of channels x 2 matrix of the channels and their
%               end nodes
%   nurbs: an array containing the nurbs description of each channel. it can be: 
%         i) an array of nurbs structures or
%         ii) number of channels x number of dimensions matrix
%             of the nodal coordinates
%   designParams:
%   options.figBounds: [xmin,xmax,ymin,ymax]
%   options.plotStyle: 1, 2
%   options.showChannels:
%   options.showPts: show end points of channels
%   options.showDesignParams: show design parameter labels
function [] = plot_channels_and_design_parameters(chan_nodes,chan_pts,nurbs,...
                                                  designParams,options)
nChans = size(chan_nodes,1);
if (numel(nurbs) ~= nChans)
    error('number of nurbs should be equal to the number of rows in chan_nodes')
end
if (nargin < 5)
    options = struct;
end
if (~isfield(options,'figBounds'))
    options.figBounds = [];
end
if (~isfield(options,'plotStyle'))
    options.plotStyle = 2;
end
if (options.plotStyle ~= 1 && options.plotStyle ~= 2)
    error('unknown plot style')
end
if (~isfield(options,'showChannels'))
    options.showChannels = true;
end
if (~isfield(options,'showPts'))
    options.showPts = false;
end
if (~isfield(options,'showDesignParams'))
    options.showDesignParams = true;
end
if (~isfield(options,'showParamLabels'))
    options.showParamLabels = true;
end
if (~isfield(options,'showParamBounds'))
    options.showParamBounds = true;
end
if ~isfield(options,'adjustFig')
    options.adjustFig = true;
end
ftsize = 20;
hold on

if (options.plotStyle == 1)
    colors = {'r','b','k','m','g','c'};
elseif (options.plotStyle == 2)
    colors = {'r'};
end
nColors = numel(colors);
specs.width = 2;
%specs.ctrlPt = 'bd';
subd = 100;

maxNCtrlPts = -inf;
for i = 1:numel(nurbs)
    specs.color = colors{rem((i-1),nColors)+1};
    specs.linestyle = '-';
    specs.width = 3;
    specs.ctrlPt = 'd';
    specs.markersize = 20;
    if (nurbs(i).number > maxNCtrlPts)
        maxNCtrlPts = nurbs(i).number;
    end
    nurbs_plot(nurbs(i),subd,specs);
    if (options.showChannels)
        midPt = nrbeval(nurbs(i),0.5);
        text(midPt(1),midPt(2),['(',num2str(i),')'],'color',specs.color,'fontsize',ftsize);
    end
end

if (options.showPts)
    for i = 1:size(chan_pts,1)
        plot(chan_pts(i,1),chan_pts(i,2),'ko','markerfacecolor','k','markersize',12)
        text(chan_pts(i,1),chan_pts(i,2),['<',num2str(i),'>'],'color','k','fontsize',ftsize);
    end
end
if (options.plotStyle == 1)
    markers = {'o','x','^','s','+','v','d'};
    lineStyles = {'-','--','-.'};
elseif (options.plotStyle == 2)
    markers = {'d'};
    lineStyles = {'-'};
    colors = {'b'};
    nColors = numel(colors);
end
nMarkers = numel(markers);
nLineStyles = numel(lineStyles);
hAlignment = {'left','left','right','right'};
nhAlignment = numel(hAlignment);
vAlignment = {'top','bottom','bottom','top'};
nvAlignment = numel(vAlignment);
xb = cell(nChans,maxNCtrlPts);
yb = cell(nChans,maxNCtrlPts);
alignInd = 0;
if (options.showDesignParams)
    for i = 1:designParams.nParams
        if (strcmpi(designParams.type{i},'CTRL_PT'))
            dim = designParams.ctrlPtDim(i);
            for j = 1:numel(designParams.channelNum{i})
                chanNum = designParams.channelNum{i}(j);
                cptNum = designParams.ctrlPtNum{i}(j);
                color = colors{rem((chanNum-1),nColors)+1};
                marker = markers{rem((j-1),nMarkers)+1};
                lineSpec = [color,marker];
                xx = nurbs(chanNum).coefs(1,cptNum)./nurbs(chanNum).coefs(4,cptNum);
                yy = nurbs(chanNum).coefs(2,cptNum)./nurbs(chanNum).coefs(4,cptNum);
                plot(xx,yy,...
                     lineSpec,'markersize',12,'markerfacecolor',color);  
                if (options.showParamLabels)
                    alignInd = alignInd + 1;           
                    hAlign = hAlignment{rem(alignInd,nhAlignment)+1};
                    vAlign = vAlignment{rem(alignInd,nvAlignment)+1};
                    text(xx,yy,...
                        ['[',num2str(i),']'],...
                        'color',color,'fontsize',ftsize,...
                        'VerticalAlignment', vAlign,...
                        'HorizontalAlignment',hAlign)

                    hAlign = hAlignment{rem(alignInd - 1,nhAlignment)+1};
                    vAlign = vAlignment{rem(alignInd - 1,nvAlignment)+1};

                    text(xx,yy,...
                        ['[',num2str(i),',',num2str(chanNum),',',num2str(cptNum),']'],...
                        'color',color,'fontsize',ftsize,...
                        'VerticalAlignment', vAlign,...
                        'HorizontalAlignment',hAlign)
                end
                if (options.showParamBounds)
                    if (dim == 1)
                        xb{chanNum,cptNum} = designParams.bounds(i,:);
                        if (isempty(yb{chanNum,cptNum}))
                            yb{chanNum,cptNum} = yy;
                        end
                    elseif (dim == 2)
                        if (isempty(xb{chanNum,cptNum}))
                            xb{chanNum,cptNum} = xx;
                        end
                        yb{chanNum,cptNum} = designParams.bounds(i,:);

                    end
                end
            end

        end
    end
end

if (options.showDesignParams && options.showParamBounds)
    if (options.plotStyle == 2)
        colors = {'k'};
        lineStyles = {'--'};
        nColors = numel(colors);
        nLineStyles = numel(lineStyles);
    end
    for i = 1:nChans
        for j = 1:maxNCtrlPts
            if (numel(xb{i,j}) == 2 && numel(yb{i,j}) == 2)   
                style = lineStyles{rem((i-1),nLineStyles)+1};
                color = colors{rem((i-1),nColors)+1}; 
                rectangle('Position',...
                          [xb{i,j}(1),yb{i,j}(1),...
                           xb{i,j}(2)-xb{i,j}(1),...
                           yb{i,j}(2)-yb{i,j}(1)],...
                          'LineStyle',style,'EdgeColor',color,'LineWidth',3)
            elseif (numel(xb{i,j}) == 2 && numel(yb{i,j}) == 1)
                style = lineStyles{rem((i-1),nLineStyles)+1};
                color = colors{rem((i-1),nColors)+1}; 
                plot(xb{i,j},[yb{i,j},yb{i,j}],[color,style],'linewidth',3)
            elseif (numel(xb{i,j}) == 1 && numel(yb{i,j}) == 2)
                style = lineStyles{rem((i-1),nLineStyles)+1};
                color = colors{rem((i-1),nColors)+1}; 
                plot([xb{i,j},xb{i,j}],yb{i,j},[color,style],'linewidth',3)
            end
        end
    end
end

if options.adjustFig
    set(gcf,'units','normalized','outerposition',[0 0 1 1])
    set(gcf,'color',[0.95,0.95,0.95])
    set(gca,'fontsize',ftsize)

    xlabel('x', 'fontsize',ftsize);
    ylabel('y', 'fontsize',ftsize);
    set(gca,'fontsize',ftsize)

    axis image

    box on
    bounds = [get(gca,'xlim'),get(gca,'ylim')];
    dels = [bounds(2)-bounds(1),bounds(4)-bounds(3)];
    maxDels = max(dels);
    newDels = [max(dels(1),0.1*maxDels),...
               max(dels(2),0.1*maxDels)];

    centers = 0.5*[bounds(1)+bounds(2),bounds(3)+bounds(4)];
    xlim([centers(1)-0.5*newDels(1),centers(1)+0.5*newDels(1)])
    ylim([centers(2)-0.5*newDels(2),centers(2)+0.5*newDels(2)])
    %axis off
    %set(gca,'units','normalized','position',[0.15,0.15,0.7,0.7])
end
end