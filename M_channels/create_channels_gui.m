function create_channels_gui
path(path,'../../NURBS/nurbs_toolbox')

%  Create and then hide the GUI as it is being constructed.
clc
close all
f = figure;

ftsize = 16;
% Construct the components.
%{
h1    = uicontrol('Style','pushbutton',...
             'String','Point 1','units','normalized','Position',[0.8,0.5,0.1,0.05],...
             'Callback',@p1button_Callback);
h2    = uicontrol('Style','pushbutton',...
             'String','Point 2','units','normalized','Position',[0.8,0.45,0.1,0.05],...
             'Callback',@p2button_Callback);
%}

% Axix bounds
uicontrol('Style','text',...
          'String','Axes bounds ',...
          'units','normalized',...
          'Position',[0.58,0.92,0.075,0.03]);
uicontrol('Style','pushbutton',...
          'String','Tight','units','normalized','Position',[0.655,0.92,0.075,0.03],...
          'Callback',@axisTightButton_Callback);      
uicontrol('Style','text',...
          'String','x1,x2,y1,y2 = ',...
          'units','normalized',...
          'Position',[0.58,0.88,0.075,0.03]);   
uicontrol('Style','edit',...
          'String','',...
          'units','normalized',...
          'Position',[0.655,0.88,0.075,0.03],...
          'BackgroundColor','w',...
          'Callback',@axisLims_Callback);  

% channel end points
endPointTitlePosWH = [0.58,0.84,0.15,0.03];
endPointNum = 1;
hEndPointNum = uicontrol('Style','text',...
                          'String',['End point = ',num2str(endPointNum)],...
                          'units','normalized',...
                          'Position',endPointTitlePosWH);
uicontrol('Style','text',...
          'String','x,y = ',...
          'units','normalized',...
          'Position',[0.58,0.8,0.075,0.03]);  
endPointPosWH = [0.655,0.8,0.075,0.03];
hEndPoint= uicontrol('Style','edit',...
                      'String','',...
                      'units','normalized',...
                      'Position',endPointPosWH,...
                      'BackgroundColor','w',...
                      'Callback',@endPoint_Callback); 
uicontrol('Style','pushbutton',...
          'String','Previous',...
          'units','normalized',...
          'Position',[0.58,0.76,0.05,0.03],...
          'Callback',@endPointPrevButton_Callback); 
uicontrol('Style','pushbutton',...
          'String','Next',...
          'units','normalized',...
          'Position',[0.63,0.76,0.05,0.03],...
          'Callback',@endPointNextButton_Callback); 
uicontrol('Style','pushbutton',...
          'String','Remove',...
          'units','normalized',...
          'Position',[0.68,0.76,0.05,0.03],...
          'Callback',@endPointRemoveButton_Callback);      
%{
uicontrol('Style','pushbutton',...
          'String','Next','units','normalized','Position',[0.85,0.68,0.08,0.03],...
          'Callback',@endPointNextButton_Callback);      
%}

% channels
channelTitlePosWH = [0.58,0.72,0.15,0.03];
channelNum = 1;
hChannelNum = uicontrol('Style','text',...
                        'String', ['Channel = ',num2str(channelNum)],...
                        'units','normalized',...
                        'Position',channelTitlePosWH);
uicontrol('Style','text',...
          'String','i,j = ',...
          'units','normalized',...
          'Position',[0.58,0.68,0.075,0.03]);  
channelContvtyPosWH = [0.655,0.68,0.075,0.03];
hChannelContvty = uicontrol('Style','edit',...
                            'String','',...
                            'units','normalized',...
                            'Position',channelContvtyPosWH,...
                            'BackgroundColor','w',...
                            'Callback',@channelContvty_Callback); 
uicontrol('Style','pushbutton',...
          'String','Add knots','units','normalized','Position',[0.58,0.64,0.05,0.03],...
          'Callback',@channelAddKnotButton_Callback);
uicontrol('Style','pushbutton',...
          'String','Add ctrl pts','units','normalized','Position',[0.63,0.64,0.05,0.03],...
          'Callback',@channelAddCtrlPtButton_Callback); 
uicontrol('Style','pushbutton',...
          'String','Update','units','normalized','Position',[0.68,0.64,0.05,0.03],...
          'Callback',@channelUpdateKnotCtrlPtButton_Callback);
uicontrol('Style','pushbutton',...
          'String','Delete knots','units','normalized','Position',[0.58,0.6,0.05,0.03],...
          'Callback',@channelDelKnotButton_Callback);
uicontrol('Style','pushbutton',...
          'String','Delete ctrl pts','units','normalized','Position',[0.63,0.6,0.05,0.03],...
          'Callback',@channelDelCtrlPtButton_Callback);       
uicontrol('Style','text',...
          'String','Diameter = ',...
          'units','normalized',...
          'Position',[0.58,0.56,0.075,0.03]);  
defaultDiam = 5e-4;
defaultDiamStr = num2str(defaultDiam);
hChannelDiams = uicontrol('Style','edit',...
                         'String',defaultDiamStr,...
                         'units','normalized',...
                         'Position',[0.655,0.56,0.075,0.03],...
                         'BackgroundColor','w',...
                         'Callback',@channelDiam_Callback);       
      
uicontrol('Style','pushbutton',...
          'String','Previous','units','normalized','Position',[0.58,0.52,0.05,0.03],...
          'Callback',@channelPrevButton_Callback); 
uicontrol('Style','pushbutton',...
          'String','Next','units','normalized','Position',[0.63,0.52,0.05,0.03],...
          'Callback',@channelNextButton_Callback); 
uicontrol('Style','pushbutton',...
          'String','Remove','units','normalized','Position',[0.68,0.52,0.05,0.03],...
          'Callback',@channelRemoveButton_Callback); 
      
% tables for knots and control points
knotColfmt = {'numeric'};
knotColwdt = {70};
knotColnames = {'knots'};
knotColedit = true;
hKnotTable = uitable('Units', 'normalized',...
        'Position', [0.75,0.02,0.07,0.92],... 
        'Data',  0,...
        'ColumnName', knotColnames,...
        'ColumnFormat', knotColfmt,...
        'ColumnWidth', knotColwdt,...
        'ColumnEditable', knotColedit,...
         'CellEditCallback', {@knotTable_Callback});

ctrlPtColfmt = repmat({'numeric'},[1,3]);
ctrlPtColwdt = repmat({75},[1,3]);
ctrlPtColnames = {'x','y','w'};
ctrlPtColedit = true(1,3);
hCtrlPtTable = uitable('Units', 'normalized',...
        'Position', [0.82,0.02,0.16,0.92],...
        'Data',  [0,0,0],...
        'ColumnName', ctrlPtColnames,...
        'ColumnFormat', ctrlPtColfmt,...
        'ColumnWidth', ctrlPtColwdt,...
        'ColumnEditable', ctrlPtColedit,...
         'CellEditCallback', {@ctrlPtTable_Callback});
 
% channel variables
channels.pts = [];
channels.contvty = [];
knots = [];
ctrlPts = [];
pt_channels = [];
channels.diams = defaultDiam;
channels.crossSection = 'square';
channels.model = 1;
channels.pt_temp = [1,27];
channels.inletEndPoint = 1;
channels.massin = 5e-4;
channels.pressureOutletEndPoint = 1;
channels.pressureOut = 0.0;
channels.heatCapacity = 3494; 
channels.viscosity = 3.405e-6;
channels.density = 1065;
% models, cross section geometry, boundary conditions and fluid properties
hModel = uicontrol('Style','popupmenu',...
                   'String',{'Mean Temperature','Constant Heat'},...
                   'units','normalized','Position',[0.58,0.48,0.075,0.03],...
                   'Callback',@modelMenu_Callback);
hCrossSection = uicontrol('Style','popupmenu',...
                          'String',{'Square','Circular'},...
                          'units','normalized','Position',[0.655,0.48,0.075,0.03],...
                          'Callback',@crossSectionMenu_Callback);
uicontrol('Style','text',...
          'String','Pt',...
          'units','normalized',...
          'Position',[0.58,0.44,0.035,0.03]);  
hTempPt = uicontrol('Style','edit',...
        'String',num2str(channels.pt_temp(1)),...
        'units','normalized',...
        'Position',[0.615,0.44,0.04,0.03],...
        'BackgroundColor','w',...
        'Callback',@channelTempPt_Callback); 
uicontrol('Style','text',...
          'String','Temp',...
          'units','normalized',...
          'Position',[0.655,0.44,0.035,0.03]);  
hPtTemp = uicontrol('Style','edit',...
        'String',num2str(channels.pt_temp(2)),...
        'units','normalized',...
        'Position',[0.69,0.44,0.04,0.03],...
        'BackgroundColor','w',...
        'Callback',@channelPtTemp_Callback);
                
uicontrol('Style','text',...
          'String','Inlet',...
          'units','normalized',...
          'Position',[0.58,0.4,0.035,0.03]);                  
hInlet = uicontrol('Style','edit',...
        'String',num2str(channels.inletEndPoint),...
        'units','normalized',...
        'Position',[0.615,0.4,0.04,0.03],...
        'BackgroundColor','w',...
        'Callback',@channelInlet_Callback); 
uicontrol('Style','text',...
          'String','Mass in',...
          'units','normalized',...
          'Position',[0.655,0.4,0.035,0.03]);                  
hMassIn = uicontrol('Style','edit',...
        'String',num2str(channels.massin),...
        'units','normalized',...
        'Position',[0.69,0.4,0.04,0.03],...
        'BackgroundColor','w',...
        'Callback',@channelMassIn_Callback); 
                
uicontrol('Style','text',...
          'String','Outlet',...
          'units','normalized',...
          'Position',[0.58,0.36,0.035,0.03]);                  
hOutlet = uicontrol('Style','edit',...
        'String',num2str(channels.pressureOutletEndPoint),...
        'units','normalized',...
        'Position',[0.615,0.36,0.04,0.03],...
        'BackgroundColor','w',...
        'Callback',@channelOutlet_Callback); 
uicontrol('Style','text',...
          'String','Pressure',...
          'units','normalized',...
          'Position',[0.655,0.36,0.035,0.03]);                  
hPressureOut = uicontrol('Style','edit',...
        'String',num2str(channels.pressureOut),...
        'units','normalized',...
        'Position',[0.69,0.36,0.04,0.03],...
        'BackgroundColor','w',...
        'Callback',@channelPressureOut_Callback);
                
                
uicontrol('Style','text',...
          'String','Heat capacity',...
          'units','normalized',...
          'Position',[0.58,0.32,0.035,0.04]);                  
hHeatCapacity = uicontrol('Style','edit',...
        'String',num2str(channels.heatCapacity),...
        'units','normalized',...
        'Position',[0.615,0.32,0.04,0.03],...
        'BackgroundColor','w',...
        'Callback',@heatCapacity_Callback); 
uicontrol('Style','text',...
          'String','Viscosity',...
          'units','normalized',...
          'Position',[0.655,0.32,0.035,0.03]);                  
hViscosity = uicontrol('Style','edit',...
        'String',num2str(channels.viscosity),...
        'units','normalized',...
        'Position',[0.69,0.32,0.04,0.03],...
        'BackgroundColor','w',...
        'Callback',@viscosity_Callback); 
uicontrol('Style','text',...
          'String','density',...
          'units','normalized',...
          'Position',[0.58,0.28,0.035,0.03]);                  
hDensity = uicontrol('Style','edit',...
        'String',num2str(channels.density),...
        'units','normalized',...
        'Position',[0.615,0.28,0.04,0.03],...
        'BackgroundColor','w',...
        'Callback',@density_Callback); 
    
% design parameters
paramNum = 1;
hParamNum = uicontrol('Style','text',...
                      'String', ['Parameter = ',num2str(paramNum)],...
                      'units','normalized',...
                      'Position',[0.58,0.22,0.075,0.03]);

hParamType = uicontrol('Style','popupmenu',...
                      'String',{'Ctrl pt','Diameter'},...
                      'units','normalized','Position',[0.655,0.22,0.075,0.03],...
                      'Callback',@paramTypeMenu_Callback);

uicontrol('Style','text',...
          'String','Channel',...
          'units','normalized',...
          'Position',[0.58,0.18,0.035,0.03]); 
hParamChannel = uicontrol('Style','edit',...
                            'String','',...
                            'units','normalized',...
                            'Position',[0.615,0.18,0.04,0.03],...
                            'BackgroundColor','w',...
                            'Callback',@paramChannel_Callback); 
uicontrol('Style','text',...
          'String','(lb,ub)',...
          'units','normalized',...
          'Position',[0.655,0.18,0.035,0.03]);
hParamBound = uicontrol('Style','edit',...
                         'String','',...
                         'units','normalized',...
                         'Position',[0.69,0.18,0.04,0.03],...
                         'BackgroundColor','w',...
                         'Callback',@paramBound_Callback); 
                     
hParamCtrlPtTitle = uicontrol('Style','text',...
                             'String','Ctrl pt',...
                             'units','normalized',...
                             'Position',[0.58,0.14,0.035,0.03]); 
hParamCtrlPt = uicontrol('Style','edit',...
                         'String','',...
                         'units','normalized',...
                         'Position',[0.615,0.14,0.04,0.03],...
                         'BackgroundColor','w',...
                         'Callback',@paramCtrlPt_Callback);                         

hParamDimTitle = uicontrol('Style','text',...
                           'String','Dim',...
                           'units','normalized',...
                           'Position',[0.655,0.14,0.035,0.03]);
hParamDim = uicontrol('Style','edit',...
                         'String','',...
                         'units','normalized',...
                         'Position',[0.69,0.14,0.04,0.03],...
                         'BackgroundColor','w',...
                         'Callback',@paramDim_Callback);
                     
          

uicontrol('Style','pushbutton',...
          'String','Previous','units','normalized','Position',[0.58,0.1,0.05,0.03],...
          'Callback',@paramPrevButton_Callback); 
uicontrol('Style','pushbutton',...
          'String','Next','units','normalized','Position',[0.63,0.1,0.05,0.03],...
          'Callback',@paramNextButton_Callback); 
uicontrol('Style','pushbutton',...
          'String','Remove','units','normalized','Position',[0.68,0.1,0.05,0.03],...
          'Callback',@paramRemoveButton_Callback); 


designParams.nParams = 0;
designParams.channelNum = cell(paramNum,1);
designParams.ctrlPtNum = cell(paramNum,1);
designParams.ctrlPtDim = nan(paramNum,1);
designParams.type{1} = 'CTRL_PT';
designParams.bounds = inf(paramNum,2);

uicontrol('Style','pushbutton',...
          'String','All diams as params','units','normalized','Position',[0.58,0.06,0.075,0.03],...
          'Callback',@allDiamParamButton_Callback); 
uicontrol('Style','text',...
          'String','(lb,ub)',...
          'units','normalized',...
          'Position',[0.655,0.06,0.035,0.03]);
diamBounds = [3.5e-4,1.5e-3];      
uicontrol('Style','edit',...
         'String',[num2str(diamBounds(1)),',',num2str(diamBounds(2))],...
         'units','normalized',...
         'Position',[0.69,0.06,0.04,0.03],...
         'BackgroundColor','w',...
         'Callback',@allDiamBound_Callback);

                     
% read output file
uicontrol('Style','pushbutton',...
          'String','Read file',...
          'units','normalized',...
          'Position',[0.58,0.02,0.075,0.03],...
          'Callback',@readFileButton_Callback);
      
% write output file     
uicontrol('Style','pushbutton',...
          'String','Write file',...
          'units','normalized',...
          'Position',[0.655,0.02,0.075,0.03],...
          'Callback',@writeFileButton_Callback);
      
% prepare axes and plot area      
hAxes = axes('Units','normalized',...
             'Position',[0.05,0.1,0.52,0.85],...
             'fontsize',ftsize,...
             'XLimMode', 'auto',...
             'YLimMode', 'auto');
xlabel('x','fontsize',ftsize)
ylabel('y','fontsize',ftsize)
hold on



% for nurbs_plot
nurbsSpecs.curve = 'r-';
nurbsSpecs.width = 2;
nurbsSpecs.ctrlPt = 'bd';
nurbPlotSubd = 100;
 

    function axisTightButton_Callback(source,eventdata) 
        set(hAxes,'XLimMode','auto','YLimMode','auto')
        axis image
    end
    function axisLims_Callback(source,eventdata,handles) 
        input = get(source,'string');
        splitStr = regexp(input,',','split');
        if (numel(splitStr) ~= 4)
           errordlg('must provide four numbers separated by comma','Invalid Input','modal')
        end
        x1 = str2double(splitStr{1});
        x2 = str2double(splitStr{2});
        y1 = str2double(splitStr{3});
        y2 = str2double(splitStr{4});
        if (isnan(x1) || isnan(x2) || isnan(y1) || isnan(y2))
           errordlg('enter four numbers','Invalid Input','modal')
           return
        else
           %xlim([x1,x2])
           %ylim([y1,y2])
           set(hAxes,'XLimMode','manual','YLimMode','manual','XLim',[x1,x2],'YLim',[y1,y2])
           %axis image
        end
    end
    
    function plot_end_points(pts)
       plot(pts(:,1),pts(:,2),'mo','markersize',10,'linewidth',2)

       text(pts(:,1),pts(:,2),...
            strtrim(cellstr(num2str((1:size(pts,1))'))),...
            'color','m',...
            'fontsize',ftsize,...
            'VerticalAlignment','Bottom')

    end
    
    function refresh_end_point_boxes()
        if (endPointNum > 0 && endPointNum <= size(channels.pts,1))
            endPointStr = [num2str(channels.pts(endPointNum,1)), ','...
                               num2str(channels.pts(endPointNum,2))];
        else
            endPointStr = '';
        end
        endPointNumStr = ['End point = ', num2str(endPointNum)];
        set(hEndPointNum,'String',endPointNumStr)
        set(hEndPoint,'String',endPointStr)      
          
    end
    
    


    function endPoint_Callback(source,eventdata,handles) 
        input = get(source,'string');
        splitStr = regexp(input,',','split');
        if (numel(splitStr) ~= 2)
           errordlg('must provide two numbers separated by comma','Invalid Input','modal')
           return
        end
        x = str2double(splitStr{1});
        y = str2double(splitStr{2});
        if (isnan(x) || isnan(y))
           errordlg('enter a pair of numbers','Invalid Input','modal')
           return
        else
           channels.pts(endPointNum,:) = [x,y];
           cla
           plot_end_points(channels.pts)
           if(isfield(channels,'nurbs'))
                plot_channels(channels.nurbs)
           end
           refresh_end_point_boxes();  
        end
    end
    


    function endPointPrevButton_Callback(source,eventdata) 
        if (endPointNum > 1)
            endPointNum = endPointNum - 1;
        end             
        refresh_end_point_boxes();     
    end

    function endPointNextButton_Callback(source,eventdata)
        if (endPointNum <= size(channels.pts,1))
            endPointNum = endPointNum + 1;
        end        
        refresh_end_point_boxes(); 
    end
    
    function endPointRemoveButton_Callback(source,eventdata) 
        if (endPointNum>=1 && endPointNum <= size(channels.pts,1))
             channels.pts(endPointNum,:) = [];
             refresh_end_point_boxes();
             cla
             plot_end_points(channels.pts)
           
             if (isfield(channels,'nurbs') && size(channels.contvty,1) > 0)
                channels.contvty(pt_channels{endPointNum},:) = [];
                channels.nurbs(pt_channels{endPointNum}) = [];   
                channels.diams(pt_channels{endPointNum}) = [];
                channels_sharing_pt(channels.contvty,endPointNum);
                plot_channels(channels.nurbs)
             end
             
        end   
    end
    

    function channels_sharing_pt(contvty,npts)
        pt_channels = cell(npts,1);
        for i = 1:size(contvty,1)
            pt_channels{contvty(i,1)}(end+1) = i;
            pt_channels{contvty(i,2)}(end+1) = i;
        end
    end
    %
    function plot_channels(nurbs)
        for i = 1:numel(nurbs)
            nurbs_plot(nurbs(i),nurbPlotSubd,nurbsSpecs);
            midPt = nrbeval(nurbs(i),0.5);
            text(midPt(1),midPt(2),...
                 ['(',num2str(i),')'],'color','r',...
                 'fontsize',ftsize,...
                 'VerticalAlignment','Bottom')
        end
    end
    %
    function refresh_channel_boxes()
        if (channelNum > 0 && channelNum <= size(channels.contvty,1))
            
            channelContvtyStr = [num2str(channels.contvty(channelNum,1)), ','...
                               num2str(channels.contvty(channelNum,2))];
            knots = channels.nurbs(channelNum).knots;
            refresh_knotTable();
            ctrlPts = channels.nurbs(channelNum).coefs;
            refresh_ctrlPtTable();
            channelDiamStr = num2str(channels.diams(channelNum));
        else
            channelContvtyStr = '';
            channelDiamStr = defaultDiamStr;
            knots = 0;
            ctrlPts = [0,0,0,0]';
        end
        channelStr = ['Channel = ', num2str(channelNum)]; 
        set(hChannelNum,'String',channelStr)
        set(hChannelContvty,'String',channelContvtyStr)
        set(hChannelDiams,'String',channelDiamStr)     
    end
    

    function channelContvty_Callback(source,eventdata,handles) 
        input = get(source,'string');
        splitStr = regexp(input,',','split');
        if (numel(splitStr) ~= 2)
           errordlg('must provide two numbers separated by comma','Invalid Input','modal')
           return
        end
        i = str2double(splitStr{1});
        j = str2double(splitStr{2});
        npts = size(channels.pts,1);
        if (isnan(i) || isnan(j))
           errordlg('enter a pair of numbers','Invalid Input','modal')
           return
        elseif (i > npts || j > npts)
           errordlg('0 <= i,j <= number of points','Invalid Input','modal')
           return 
        else
           channels.contvty(channelNum,:) = [i,j];
           knots = [0,0,1,1];
           ctrlPts = [channels.pts(i,1),channels.pts(j,1);
                      channels.pts(i,2),channels.pts(j,2);
                      0,0;
                      1,1];
           channels.nurbs(channelNum) = nrbmak(ctrlPts,knots);       
    
          
          
           channels_sharing_pt(channels.contvty,endPointNum);
           cla           
           plot_channels(channels.nurbs)
           plot_end_points(channels.pts)
           
         
           channelDiam_Callback(hChannelDiams,[]);
           refresh_channel_boxes();
        end
    end
    function channelDiam_Callback(source,eventdata)
        input = str2num(get(source,'string'));
        if (numel(input) ~= 1 || isnan(input))
           errordlg('must provide one number','Invalid Input for diameter','modal')
           return
        end
        if (channelNum > 0)
            channels.diams(channelNum) = input;
        end       
    end
    function channelPrevButton_Callback(source,eventdata) 
        if (channelNum > 1)
            channelNum = channelNum - 1;
        end           
        refresh_channel_boxes();     
    end
  
    function channelNextButton_Callback(source,eventdata)
        if (channelNum <= size(channels.contvty,1))
            channelNum = channelNum + 1;
        end
        refresh_channel_boxes(); 
    end

    function channelRemoveButton_Callback(source,eventdata)
        if (channelNum > 0 && channelNum <= size(channels.contvty,1))
             channels.contvty(channelNum,:) = [];
             channels.nurbs(channelNum) = [];
             channels.diams(channelNum) = [];
             refresh_channel_boxes();
             cla
             plot_end_points(channels.pts);
             plot_channels(channels.nurbs);
             channels_sharing_pt(channels.contvty,endPointNum);
             
        end   
    end

    function refresh_knotTable()
        set(hKnotTable,'Data',knots')     
    end

    function refresh_ctrlPtTable()    
        set(hCtrlPtTable,'Data',ctrlPts([1,2,4],:)')
    end

    function knotTable_Callback(source, eventdata)
        knots = get(source,'Data')';

    end

    function ctrlPtTable_Callback(source, eventdata)
        tableData = get(source,'Data');
        tableData(:,1) = tableData(:,1).*tableData(:,3);
        tableData(:,2) = tableData(:,2).*tableData(:,3);
        clear ctrlPts;
        ctrlPts([1,2,4],:) = tableData';
         
    end
    
    function channelAddKnotButton_Callback(source,eventdata,handles)
        tableData = get(hKnotTable,'Data');
        tableData(end+1) = tableData(end);   
        set(hKnotTable,'Data',tableData)
    end
    
    function channelDelKnotButton_Callback(source,eventdata,handles)
        tableData = get(hKnotTable,'Data');
        tableData(end) = []; 
        set(hKnotTable,'Data',tableData)
    end

    function channelAddCtrlPtButton_Callback(source,eventdata,handles)
        tableData = get(hCtrlPtTable,'Data');
        tableData(end+1,:) = tableData(end,:);
        tableData(end-1,:) = [0,0,1];
        set(hCtrlPtTable,'Data',tableData)
    end
    
    function channelDelCtrlPtButton_Callback(source,eventdata,handles)
        tableData = get(hCtrlPtTable,'Data');
        tableData(end,:) = [];
        set(hCtrlPtTable,'Data',tableData)
    end
    function channelUpdateKnotCtrlPtButton_Callback(source,eventdata)
        knots = get(hKnotTable,'Data')';
        clear ctrlPts;
        ctrlPts([1,2,4],:) = get(hCtrlPtTable,'Data')';
        channels.nurbs(channelNum) = nrbmak(ctrlPts,knots);  
        cla
        plot_channels(channels.nurbs)
        plot_end_points(channels.pts)
    end
    
    function refresh_properties()
        set(hModel,'value',channels.model)
        if (strcmpi(channels.crossSection,'square'))
            set(hCrossSection,'value',1)
        elseif (strcmpi(channels.crossSection,'circular'))
            set(hCrossSection,'value',2)
        else
            error('unknown cross section')
        end
        if (isempty(channels.pt_temp))
            set(hTempPt,'String','NA')
            set(hPtTemp,'String','NA')
        else
            set(hTempPt,'String',num2str(channels.pt_temp(1)))
            set(hPtTemp,'String',num2str(channels.pt_temp(2)))
        end
        set(hInlet,'String',num2str(channels.inletEndPoint))
        set(hMassIn,'String',num2str(channels.massin))
        set(hOutlet,'String',num2str(channels.pressureOutletEndPoint))
        set(hPressureOut,'String',num2str(channels.pressureOut))
        set(hHeatCapacity,'String',num2str(channels.heatCapacity))
        set(hViscosity,'String',num2str(channels.viscosity))
         set(hDensity,'String',num2str(channels.density))
    end
    function modelMenu_Callback(source,eventdata)
        str = get(source, 'String');
        val = get(source,'Value');
        switch str{val}
            case 'Mean Temperature'
                channels.model = 1;
            case 'Constant Heat'
                channels.model = 2;
                set(hTempPt,'String','NA')
                set(hPtTemp,'String','NA')
            otherwise
                error('unknown model')
        end
    end

    function crossSectionMenu_Callback(source,eventdata)
        str = get(source, 'String');
        val = get(source,'Value');
        switch str{val}
            case 'Square'
                channels.crossSection = 'square';
            case 'Circular'
                channels.crossSection = 'circular';
            case 'Rectangular'
                channels.crossSection = 'rectangular';
            otherwise
                error('unknown cross section')
        end   
    end

    function channelTempPt_Callback(source,eventdata)
        input = str2num(get(source,'string'));
        if (numel(input) ~= 1 || isnan(input))
           errordlg('must provide one number','Invalid Input for pt','modal')
        else
            channels.pt_temp(1) = input;
        end
    end


    function channelPtTemp_Callback(source,eventdata)
        input = str2num(get(source,'string'));
        if (numel(input) ~= 1 || isnan(input))
           errordlg('must provide one number','Invalid Input for temp','modal')
        else
            channels.pt_temp(2) = input;
        end
    end

    function channelInlet_Callback(source,eventdata)
        input = str2num(get(source,'string'));
        if (numel(input) ~= 1 || isnan(input))
           errordlg('must provide one number','Invalid Input for inlet','modal')
        else
            channels.inletEndPoint = input;
        end
    end

    function channelMassIn_Callback(source,eventdata)
        input = str2num(get(source,'string'));
        if (numel(input) ~= 1 || isnan(input))
           errordlg('must provide one number','Invalid Input for mass in','modal')
        else
            channels.massin = input;
        end
    end

    function channelOutlet_Callback(source,eventdata)
        input = str2num(get(source,'string'));
        if (numel(input) ~= 1 || isnan(input))
           errordlg('must provide one number','Invalid Input for outlet','modal')
        else
            channels.pressureOutletEndPoint = input;
        end
    end

    function channelPressureOut_Callback(source,eventdata)
        input = str2num(get(source,'string'));
        if (numel(input) ~= 1 || isnan(input))
           errordlg('must provide one number','Invalid Input for pressure','modal')
        else
            channels.pressureOut = input;
        end
    end


    function heatCapacity_Callback(source,eventdata)
        input = str2num(get(source,'string'));
        if (numel(input) ~= 1 || isnan(input))
           errordlg('must provide one number','Invalid Input for heat capacity','modal')
        else
            channels.heatCapacity = input;
        end
    end


    function viscosity_Callback(source,eventdata)
        input = str2num(get(source,'string'));
        if (numel(input) ~= 1 || isnan(input))
           errordlg('must provide one number','Invalid Input for viscosity','modal')
        else
            channels.viscosity = input;
        end
    end

    function density_Callback(source,eventdata)
        input = str2num(get(source,'string'));
        if (numel(input) ~= 1 || isnan(input))
           errordlg('must provide one number','Invalid Input for density','modal')
        else
            channels.density = input;
        end
    end
    
    function refresh_param_boxes()
        if (designParams.nParams == 0)
            return
        end
        set(hParamNum,'String',['Parameter = ',num2str(paramNum)])
        if(isempty(designParams.channelNum{paramNum}))
            set(hParamChannel,'String','')
        else
            set(hParamChannel,'String',num2str(designParams.channelNum{paramNum}(1)))
        end
        if (any(isinf(designParams.bounds(paramNum,:))))
            set(hParamBound,'String','')
        else
            set(hParamBound,'String',[num2str(designParams.bounds(paramNum,1)),',',...
                                      num2str(designParams.bounds(paramNum,2))]);
        end
        if (strcmpi(designParams.type{paramNum},'DIAM'))
            set(hParamType,'value',2)
            set(hParamCtrlPtTitle,'String','')
            set(hParamDimTitle,'String','')
            set(hParamCtrlPt,'String','')
            set(hParamDim,'String','')
        elseif (strcmpi(designParams.type{paramNum},'CTRL_PT'))
            set(hParamType,'value',1)
            set(hParamCtrlPtTitle,'String','Ctrl pt')
            set(hParamDimTitle,'String','Dim')
            if(isempty(designParams.ctrlPtNum{paramNum}))
                set(hParamCtrlPt,'String','')
            else
                set(hParamCtrlPt,'String',num2str(designParams.ctrlPtNum{paramNum}(1)))
            end
            if(isnan(designParams.ctrlPtDim(paramNum)))
                set(hParamDim,'String','')
            else
                set(hParamDim,'String',num2str(designParams.ctrlPtDim(paramNum)))
            end
        else
            error('unknown design parameter type')
        end
    end
    
    function plot_channels_parameters()
        if (designParams.nParams == 0)
            cla
            plot_end_points(channels.pts)
            plot_channels(channels.nurbs)
            return
        end
        if ((~isempty(designParams.channelNum{paramNum}) ...
            && all(~isinf(designParams.bounds(paramNum,:)))) ...
            && (strcmpi(designParams.type{paramNum},'DIAM') ...
                || (strcmpi(designParams.type{paramNum},'CTRL_PT') ...
                && ~isempty(designParams.ctrlPtNum{paramNum}) ...
                && ~isnan(designParams.ctrlPtDim(paramNum)))))
            cla
            plot_end_points(channels.pts)
            options.showDesignParams = true;
            options.showParamBounds = true;
            options.showParamLabels = true;
            options.adjustFig = false; % must be false to avoid screwing with the buttons 
            plot_channels_and_design_parameters(channels.contvty,...
                                                channels.pts,...
                                                channels.nurbs,...
                                                designParams, ...
                                                options)
        end
    end

    function number_of_design_parameters()
        %designParams.nParams = min([nnz(~cellfun('isempty',designParams.channelNum)),...
        %                            nnz(~cellfun('isempty',designParams.ctrlPtNum)),...
        %                            nnz(~isnan(designParams.ctrlPtDim)),...
        %                            nnz(~any(isinf(designParams.bounds),2))]);
        designParams.nParams = min([nnz(~cellfun('isempty',designParams.channelNum)),...
                                    nnz(~any(isinf(designParams.bounds),2))]);
        
    end
    function paramTypeMenu_Callback(source,eventdata)
        str = get(source, 'String');
        val = get(source,'Value');
        switch str{val}
            case 'Ctrl pt'
                set(hParamCtrlPtTitle,'String','Ctrl pt')
                set(hParamDimTitle,'String','Dim')
                designParams.type{paramNum} = 'CTRL_PT';
            case 'Diameter'
                set(hParamCtrlPtTitle,'String','')
                set(hParamDimTitle,'String','')
                designParams.type{paramNum} = 'DIAM';
            otherwise
                error('unknown design parameter type')
        end
        
    end

    function paramChannel_Callback(source,eventdata)
        input = str2num(get(source,'string'));
        if (numel(input) ~= 1 || isnan(input))
           errordlg('must provide one number for channel number','Invalid Input for channel number','modal')
           return
        end
        if (paramNum > 0)
            designParams.channelNum{paramNum} = input;
            number_of_design_parameters();
            plot_channels_parameters();
        end    
    end
    
    function paramBound_Callback(source,eventdata)
        input = get(source,'string');
        splitStr = regexp(input,',','split');
        if (numel(splitStr) ~= 2)
           errordlg('must provide two numbers separated by comma for bound','Invalid Input for bound','modal')
           return
        end
        lb = str2double(splitStr{1});
        ub = str2double(splitStr{2});
        if (isnan(lb) || isnan(ub))
           errordlg('enter a pair of numbers for param bound','Invalid Input','modal')
        elseif (paramNum > 0)
           designParams.bounds(paramNum,1) = lb;
           designParams.bounds(paramNum,2) = ub;
           number_of_design_parameters();
           plot_channels_parameters();
        end
    end

    function paramCtrlPt_Callback(source,eventdata)
        if (strcmpi(designParams.type{paramNum},'CTRL_PT'))
            input = str2num(get(source,'string'));
            if (numel(input) ~= 1 || isnan(input))
               errordlg('must provide one number for param ctrl pt num','Invalid Input for param ctrl pt num','modal')
               return
            end
            if (paramNum > 0)
                designParams.ctrlPtNum{paramNum} = input;
                number_of_design_parameters();
                plot_channels_parameters();
            end  
        end
    end
    
    function paramDim_Callback(source,eventdata)
        if (strcmpi(designParams.type{paramNum},'CTRL_PT'))
            input = str2num(get(source,'string'));
            if (numel(input) ~= 1 || isnan(input))
               errordlg('must provide one number for param ctrl pt dim','Invalid Input for param ctrl pt dim','modal')
               return
            end
            if (paramNum > 0)
                designParams.ctrlPtDim(paramNum) = input;
                number_of_design_parameters();
                plot_channels_parameters();
            end  
        end
    end

    
    function paramPrevButton_Callback(source,eventdata)
        if (paramNum > 1)
            paramNum = paramNum - 1;
        end         
       
        refresh_param_boxes();
        plot_channels_parameters();
    end
    
    function paramNextButton_Callback(source,eventdata)
        if (paramNum < designParams.nParams)
            paramNum = paramNum + 1;
        elseif (paramNum == designParams.nParams)
            % check input
            moveOn = false;
            if (strcmpi(designParams.type{paramNum},'CTRL_PT'))
                if (~isempty(designParams.ctrlPtNum{paramNum}) ...
                    && ~isnan(designParams.ctrlPtDim(paramNum)) ...
                    && ~isempty(designParams.channelNum{paramNum}) ...
                    && ~any(isinf(designParams.bounds(paramNum,:))))
                    moveOn = true;
                end
            elseif (strcmpi(designParams.type{paramNum},'DIAM'))
                if (~isempty(designParams.channelNum{paramNum}) ...
                    && ~any(isinf(designParams.bounds(paramNum,:))))
                    moveOn = true;
                end
            end
            if (moveOn)
                paramNum = paramNum + 1;
                designParams.channelNum{paramNum} = [];
                designParams.ctrlPtNum{paramNum} = [];
                designParams.ctrlPtDim(paramNum) = nan;
                designParams.type{paramNum} = 'CTRL_PT';
                designParams.bounds(paramNum,:) = inf;
            end
        else        
            paramChannel_Callback(hParamChannel);
            paramBound_Callback(hParamBound);
            paramCtrlPt_Callback(hParamCtrlPt);
            paramDim_Callback(hParamDim);
        end
        refresh_param_boxes(); 
        plot_channels_parameters(); 
    end
    

    function paramRemoveButton_Callback(source,eventdata)
        if (paramNum > 0 && paramNum <= designParams.nParams)
             designParams.channelNum(paramNum) = [];
             designParams.ctrlPtNum(paramNum) = [];
             designParams.ctrlPtDim(paramNum) = [];
             designParams.type(paramNum) = [];
             designParams.bounds(paramNum,:) = [];  
             designParams.nParams = designParams.nParams - 1;
             refresh_param_boxes();
             plot_channels_parameters();
             %plot_end_points(channels.pts);
             %plot_channels(channels.nurbs);        
        end 
    end

    function allDiamParamButton_Callback(source,eventdata)
        choice = questdlg('Generate design parameter for every channel diameter?',...
                          'Diameters as design parameters','OK','Cancel','Cancel');
        switch choice
            case 'OK'
                if (any(isnan(diamBounds)))
                     errordlg('must provide two numbers separated by comma for diam bound','Invalid Input for diam bound','modal')
                else
                    nChannels = numel(channels.nurbs);
                    is = designParams.nParams + 1;
                    ie = designParams.nParams + nChannels;
                    [designParams.type{is:ie}] = deal('DIAM');
                   
                    designParams.channelNum(is:ie) = cell(1,1);
                    for i = 1:nChannels
                        %designParams.type{i + designParams.nParams} = 'DIAM';
                        designParams.channelNum{i + designParams.nParams} = i;
                    end
                    designParams.bounds(is:ie,:) = repmat(diamBounds,nChannels,1);
                    designParams.nParams = ie;
                    paramNum = designParams.nParams;
                    refresh_param_boxes();
                end
            case 'Cancel'
        end
    end

    function allDiamBound_Callback(source,eventdata)
        input = get(source,'string');
        splitStr = regexp(input,',','split');
        if (numel(splitStr) ~= 2)
           errordlg('must provide two numbers separated by comma for diam bound','Invalid Input for diam bound','modal')
           return
        end
        lb = str2double(splitStr{1});
        ub = str2double(splitStr{2});
        if (isnan(lb) || isnan(ub))
           errordlg('enter a pair of numbers for diam bound','Invalid Input','modal')
        elseif (paramNum > 0)
           diamBounds = [lb,ub];
        end
    end
    function readFileButton_Callback(source,eventdata)
        [fileName,pathName] = uigetfile('*','Open','../ChannelFiles/');
        readOptions.guimode = true;
        [channels,designParams] = preprocess_channels([pathName,fileName],readOptions);
        if (designParams.nParams == 0)
            designParams.type{1} = 'CTRL_PT';
        end
        endPointNum = size(channels.pts,1);
        refresh_end_point_boxes()
        channelNum = numel(channels.nurbs);      
        refresh_channel_boxes()
        channels_sharing_pt(channels.contvty,endPointNum);
        knots = channels.nurbs(channelNum).knots;
        refresh_knotTable();
        ctrlPts = channels.nurbs(channelNum).coefs;
        refresh_ctrlPtTable();
        refresh_properties();
        paramNum = designParams.nParams;
        if (designParams.nParams == 0)
            paramNum = 1;
        end
        refresh_param_boxes();
        figure(f)
        plot_channels_parameters();
       
    end
    function writeFileButton_Callback(source,eventdata)
        if (isfield(channels,'nurbs'))
            [filename,pathname] = uiputfile('*.channel','Save as','../ChannelFiles/');
            write_channel_file([pathname,filename],...
                                channels,...
                                designParams,...
                                'w',[]);
        end
    end
end