function CROP_tool(varargin);
%function CROP_tool(varargin);
% Crop Tool to be used with imagescn. 
% Usage: CROP_tool;
%

if isempty(varargin) 
   Action = 'New';
else
   Action = varargin{1};  
end

switch Action
case 'New'
    Create_New_Button;

case 'Activate_Crop'
    Activate_Crop;
    
case 'Deactivate_Crop'
    Deactivate_Crop(varargin{2:end});

case 'CROP_Draw'
    CROP_Draw(gcbo);
    
case 'CROP_Save'
    CROP_Save(gcbo);
    
case 'CROP_New_Figure'
    CROP_New_Figure;
    
case 'Menu_Crop'
    Menu_Crop;
    
case 'Close_Parent_Figure'
    Close_Parent_Figure;
    
otherwise
    disp(['Unimplemented Functionality: ', Action]);
   
end;
      
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function Create_New_Button
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

fig = gcf;

% Find handle for current image toolbar and menubar
hToolbar = findall(fig, 'type', 'uitoolbar', 'Tag','FigureToolBar' );
hToolMenu = findall(fig, 'Label', '&Tools');

if ~isempty(hToolbar) & isempty(findobj(hToolbar, 'Tag', 'figCropTool'))
	hToolbar_Children = get(hToolbar, 'Children');

	% The default button size is 15 x 16 x 3. Create Button Image
   button_size_x= 16;
   button_image = NaN* zeros(15,button_size_x);
   
   f=[...
    26    40    41    55    56    70    71    80    81    82    83    84    85,...
    86    87    88    89    90    95    96    97    98    99   100   101   102,...
   103   104   110   111   125   126   130   131   140   141   145   146   152,...
   153   154   155   156   158   159   160   161   166   167   168   169   170,...
   171   173   174   175   176   185   186   200   201   215   216   230];
   button_image(f) = 0;
   button_image = repmat(button_image, [1,1,3]);

   buttontags = {'figWindowLevel', 'figPanZoom', 'figROITool', 'figViewImages', 'figPointTool', 'figRotateTool', 'figProfileTool','figCropTool'};
   separator = 'off';
   
   hbuttons = [];
   for i = 1:length(buttontags)
       hbuttons = [hbuttons, findobj(hToolbar_Children, 'Tag', buttontags{i})];
   end;
   if isempty(hbuttons)
       separator = 'on';
   end;
   
   hNewButton = uitoggletool(hToolbar);
   set(hNewButton, 'Cdata', button_image, ...
      'OnCallback', 'CROP_tool(''Activate_Crop'')',...
      'OffCallback', 'CROP_tool(''Deactivate_Crop'')',...
      'Tag', 'figCropTool', ...
      'TooltipString', 'Crop figure',...
      'UserData', [], ...
      'Enable', 'on');   
end;

% If the menubar exists, create menu item
if ~isempty(hToolMenu) & isempty(findobj(hToolMenu, 'Tag', 'menuCrop'))
  hWindowLevelMenu = findobj(hToolMenu, 'Tag', 'menuWindowLevel');
  hPanZoomMenu     = findobj(hToolMenu, 'Tag', 'menuCrop');
  hROIToolMenu     = findobj(hToolMenu, 'Tag', 'menuROITool');
  hViewImageMenu   = findobj(hToolMenu, 'Tag', 'menuViewImages');
  hPointToolMenu   = findobj(hToolMenu, 'Tag', 'menuPointTool');
  hRotateToolMenu  = findobj(hToolMenu, 'Tag', 'menuRotateTool');
  hProfileToolMenu = findobj(hToolMenu, 'Tag', 'menuProfileTool');
  hCropToolMenu    = findobj(hToolMenu, 'Tag', 'menuCropTool');
	
  position = 10;
  separator = 'On';
  hMenus = [ hWindowLevelMenu, hROIToolMenu, hViewImageMenu, hPointToolMenu, hRotateToolMenu,hProfileToolMenu,hCropToolMenu ];

  if length(hMenus>0) 
	  position = position + length(hMenus);
	  separator = 'Off';
  end;
  
  hNewMenu = uimenu(hToolMenu,'Position', position);
  set(hNewMenu, 'Tag', 'menuCrop','Label',...
      'Crop',...
      'CallBack', 'CROP_tool(''Menu_Crop'')',...
      'Separator', separator,...
      'UserData', hNewButton...
  ); 
  
end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function Activate_Crop(varargin);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin ==0
    set(0, 'ShowHiddenHandles', 'On');
    hNewButton = gcbo;
    set(findobj('Tag', 'menuCrop'),'checked', 'on');
else
    hNewButton = varargin{1};
end;

% allows for calls from buttons other than those in toolbar
fig = get(hNewButton, 'Parent');
if ~strcmp(get(fig, 'Type'), 'figure'),
    fig = get(fig, 'Parent');
end


% Deactivate zoom and rotate buttons
hToolbar = findall(fig, 'type', 'uitoolbar');
hToolbar = findobj(hToolbar, 'Tag', 'FigureToolBar');
if ~isempty(hToolbar)
	hToolbar_Children = get(hToolbar, 'Children');
	
	% disable MATLAB's own tools
	Rot3D = findobj(hToolbar_Children,'Tag', 'figToolRotate3D');
	ZoomO = findobj(hToolbar_Children,'Tag', 'figToolZoomOut');
	ZoomI = findobj(hToolbar_Children,'Tag', 'figToolZoomIn');

	% try to disable other tools buttons - if they exist
	WL = findobj(hToolbar_Children, 'Tag', 'figWindowLevel');
	PZ = findobj(hToolbar_Children,'Tag', 'figPanZoom');
	RT = findobj(hToolbar_Children,'Tag', 'figROITool');
	MV = findobj(hToolbar_Children,'Tag', 'figViewImages');
	PM = findobj(hToolbar_Children,'Tag', 'figPointTool');
	RotT = findobj(hToolbar_Children,'Tag', 'figRotateTool');
	Prof = findobj(hToolbar_Children, 'Tag', 'figProfileTool');	
    CROP = findobj(hToolbar_Children,'Tag', 'figCropTool');

	old_ToolHandles  =     [Rot3D, ZoomO, ZoomI,WL,PZ, RT,MV,PM, RotT,Prof];
	old_ToolEnables  = get([Rot3D, ZoomO, ZoomI,WL,PZ, RT,MV,PM, RotT,Prof], 'Enable');
	old_ToolStates   = get([Rot3D, ZoomO, ZoomI,WL,PZ, RT,MV,PM, RotT,Prof], 'State');
	
	for i = 1:length(old_ToolHandles)
		if strcmp(old_ToolStates(i) , 'on')			
			set(old_ToolHandles(i), 'State', 'Off');
		end;
		set(old_ToolHandles(i), 'Enable', 'Off');
	end;
        %LFG
        %enable save_prefs tool button
        SP = findobj(hToolbar_Children, 'Tag', 'figSavePrefsTool');
        set(SP,'Enable','On');
end;


% Start CROP GUI
fig2_old = findobj('Tag', 'CROP_figure');
% close the old WL figure to avoid conflicts
if ~isempty(fig2_old) close(fig2_old);end;

% open new figure
fig2_file = 'CROP_tool_figure.fig';
fig2 = openfig(fig2_file,'reuse');
optional_uicontrols = { ...
    'crop_variable_name_edit', 'String'; ...
                   };
set(SP,'Userdata',{fig2, fig2_file, optional_uicontrols});

% Generate a structure of handles to pass to callbacks, and store it. 
handlesCROP = guihandles(fig2);
guidata(fig2,handlesCROP);

close_str = [ 'hNewButton = findobj(''Tag'', ''figCropTool'');' ...
        ' if strcmp(get(hNewButton, ''Type''), ''uitoggletool''),'....
        ' set(hNewButton, ''State'', ''off'' );' ...
        ' else,  ' ...
        ' CROP_tool(''Deactivate_Crop'',hNewButton);'...
        ' set(hNewButton, ''Value'', 0);',...
        ' end;' ];

set(fig2, 'Name', 'CROP Tool',...
    'closerequestfcn', close_str);

% Record and store previous WBDF etc to restore state after CROP is done. 
% old_WBDF = get(fig, 'WindowButtonDownFcn');
% old_WBMF = get(fig, 'WindowButtonMotionFcn');
% old_WBUF = get(fig, 'WindowButtonUpFcn');
old_UserData = get(fig, 'UserData');
old_CRF = get(fig, 'Closerequestfcn');

% Store initial state of all axes in current figure for reset
h_all_axes = findobj(fig,'Type','Axes');
for i = 1:length(h_all_axes)
    all_xlims(i,:) = get(h_all_axes(i),'Xlim');
    all_ylims(i,:) = get(h_all_axes(i),'Ylim');
end;

h_axes = h_all_axes(end);
set(fig, 'CurrentAxes', h_axes);

% Draw faster and without flashes
set(fig, 'Closerequestfcn', [ old_CRF , ',CROP_tool(''Close_Parent_Figure'')']);
set(fig, 'Renderer', 'zbuffer');
set(0, 'ShowHiddenHandles', 'On', 'CurrentFigure', fig);
set(gca,'Drawmode', 'Fast');

% store the figure's old infor within the fig's own userdata
set(fig, 'UserData', {fig2,  old_UserData,...
         old_CRF, old_ToolEnables,old_ToolHandles, old_ToolStates});

% store all relevant info for faster use during calls
set(handlesCROP.Title_text, 'UserData', {fig, fig2, h_all_axes, all_xlims, all_ylims, h_axes });

% CROP_Draw(handlesCROP.crop_draw_pb);

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function Deactivate_Crop(varargin);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin ==0
    set(0, 'ShowHiddenHandles', 'On');    
    hNewButton = gcbo;
    set(findobj('Tag', 'menuCrop'),'checked', 'Off');
else
    hNewButton = varargin{1};
end;
    
% Reactivate other buttons
fig = get(hNewButton, 'Parent');
if ~strcmp(get(fig, 'Type'), 'figure'),
    fig = get(fig, 'Parent');
end

hToolbar = findall(fig, 'type', 'uitoolbar');
if ~isempty(hToolbar)
    hToolbar_Children = get(hToolbar, 'Children');
    set(findobj(hToolbar_Children,'Tag', 'figToolRotate3D'),'Enable', 'On');
    set(findobj(hToolbar_Children,'Tag', 'figToolZoomOut'),'Enable', 'On');
    set(findobj(hToolbar_Children,'Tag', 'figToolZoomIn'),'Enable', 'On');

end;

myzoom('off');

% Restore old BDFs
old_info= get(fig,'UserData');

% Restore old Pointer and UserData
set(fig, 'UserData', old_info{2});
set(fig, 'CloseRequestFcn', old_info{3});
old_ToolEnables  = old_info{4};
old_ToolHandles = old_info{5};
old_ToolStates = old_info{6};

fig2 = old_info{1};

hCropbox=getappdata(findobj(fig2,'tag','crop_draw_pb'),'hCropbox');
if ~isempty(hCropbox) & ishandle(hCropbox)
    delete(hCropbox)
end

try
	set(fig2, 'CloseRequestFcn', 'closereq');
	close(fig2); 
catch
	delete(fig2);
end;    

for i = 1:length(old_ToolHandles)
	try
		set(old_ToolHandles(i), 'Enable', old_ToolEnables{i}, 'State', old_ToolStates{i});	catch
	end;
end;
%LFG
%disable save_prefs tool button
SP = findobj(hToolbar_Children, 'Tag', 'figSavePrefsTool');
set(SP,'Enable','Off');

set(0, 'ShowHiddenHandles', 'Off');



%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function CROP_Draw(hObject);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
hCropbox=getappdata(hObject,'hCropbox');
if ~isempty(hCropbox)
    delete(hCropbox(ishandle(hCropbox)))
end
k = waitforbuttonpress;
point1 = get(gca,'CurrentPoint');    % button down detected
finalRect = rbbox;                   % return figure units
point2 = get(gca,'CurrentPoint');    % button up detected
point1 = point1(1,1:2);              % extract x and y
point2 = point2(1,1:2);
p1 = min(point1,point2);             % calculate locations
offset = abs(point1-point2);         % and dimensions
x = [p1(1) p1(1)+offset(1) p1(1)+offset(1) p1(1) p1(1)];
y = [p1(2) p1(2) p1(2)+offset(2) p1(2)+offset(2) p1(2)];
handlesCROP=guidata(hObject);
userdata=get(handlesCROP.Title_text, 'UserData');
h_all_axes=userdata{3};

for i=1:length(h_all_axes)
    axes(h_all_axes(i));
    holdstate=ishold;
    hold on
    hCropbox(i)=plot(x,y,'r','linewidth',1);% draw box around selected region
    if ~holdstate
        hold off
    end
end
setappdata(hObject,'hCropbox',hCropbox);
set(hObject,'string','redraw');

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function CROP_Save(hObject);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
handlesCROP=guidata(hObject);
userdata=get(handlesCROP.Title_text, 'UserData');
h_all_axes=userdata{3};

hCropbox=getappdata(handlesCROP.crop_draw_pb,'hCropbox');
if isempty(hCropbox) | any(~ishandle(hCropbox));
    return;
else
    xdata=get(hCropbox(1),'xdata');
    ydata=get(hCropbox(1),'ydata');
    xlim=round(unique(xdata));
    ylim=round(unique(ydata));
end


if isappdata(h_all_axes(1),'ImageData')
    imagedata=getappdata(h_all_axes(1),'ImageData');
    [nframes]=size(imagedata,3);
    [xlim,ylim]=CROP_limitcheck(xlim,ylim,imagedata);
    cropdata=zeros(diff(ylim)+1,diff(xlim)+1,nframes,length(h_all_axes));
else
    h_image=findobj(h_all_axes(1),'type','image');
    imagedata=get(h_image,'cdata');
    [xlim,ylim]=CROP_limitcheck(xlim,ylim,imagedata);
    cropdata=zeros(diff(ylim)+1,diff(xlim)+1,length(h_all_axes));
end

for i=1:length(h_all_axes)
    if isappdata(h_all_axes(i),'ImageData')
        imagedata=getappdata(h_all_axes(i),'ImageData');
        cropdata(:,:,:,i)=imagedata(ylim(1):ylim(2),xlim(1):xlim(2),:);
    else
        h_image=findobj(h_all_axes(i),'type','image');
	    imagedata=get(h_image,'cdata');
        cropdata(:,:,i)=imagedata(ylim(1):ylim(2),xlim(1):xlim(2),:);
    end
end

% save crop data to workspace
assignin('base',get(handlesCROP.crop_variable_name_edit,'string'),cropdata);

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function [xlim,ylim]=CROP_limitcheck(xlim,ylim,imagedata);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
[nrows,ncols,nframes]=size(imagedata);
if xlim(1)<1;     xlim(1)=1;     end
if xlim(2)>ncols; xlim(2)=ncols; end
if ylim(1)<1;     ylim(1)=1;     end
if ylim(2)>nrows; ylim(2)=nrows; end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function CROP_New_Figure;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function Menu_Crop;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

hNewMenu = gcbo;
checked=  umtoggle(hNewMenu);
hNewButton = get(hNewMenu, 'userdata');

if ~checked
    % turn off button
    %Deactivate_Crop(hNewButton);
    set(hNewMenu, 'Checked', 'off');
    set(hNewButton, 'State', 'off' );
else
    %Activate_Crop(hNewButton);
    set(hNewMenu, 'Checked', 'on');
    set(hNewButton, 'State', 'on' );
end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function Close_Parent_Figure;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to make sure that if parent figure is closed, 
% the ROI info and ROI Tool are closed too.
set(findobj('Tag', 'CROP_figure'), 'Closerequestfcn', 'closereq');
try 
    close(findobj('Tag','CROP_figure'));
end;

