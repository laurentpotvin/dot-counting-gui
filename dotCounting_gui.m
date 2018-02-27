function varargout = dotCounting_gui(varargin)
% DOTCOUNTING_GUI MATLAB code for dotCounting_gui.fig
%      DOTCOUNTING_GUI, by itself, creates a new DOTCOUNTING_GUI or raises the existing
%      singleton*.
%
%      H = DOTCOUNTING_GUI returns the handle to a new DOTCOUNTING_GUI or the handle to
%      the existing singleton*.
%
%      DOTCOUNTING_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DOTCOUNTING_GUI.M with the given input arguments.
%
%      DOTCOUNTING_GUI('Property','Value',...) creates a new DOTCOUNTING_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before dotCounting_gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to dotCounting_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help dotCounting_gui

% Last Modified by GUIDE v2.5 24-Dec-2017 11:59:14

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @dotCounting_gui_OpeningFcn, ...
    'gui_OutputFcn',  @dotCounting_gui_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before dotCounting_gui is made visible.
function dotCounting_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to dotCounting_gui (see VARARGIN)

% Choose default command line output for dotCounting_gui
handles.output = hObject;



% Initializing variables
handles.open_stack = 0;
% handles.num_channels = 3;
handles.stack.autothreshold = 1;
handles.list_of_open_objects = []; % unused
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes dotCounting_gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);



% --- Outputs from this function are returned to the command line.
function varargout = dotCounting_gui_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --------------------------------------------------------------------
function menu_File_Callback(hObject, eventdata, handles)
% hObject    handle to menu_File (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --------------------------------------------------------------------
function menu_cells_Callback(hObject, eventdata, handles)
% hObject    handle to menu_cells (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --------------------------------------------------------------------
function menu_window_Callback(hObject, eventdata, handles)
% hObject    handle to menu_window (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --------------------------------------------------------------------
function menu_stackViewer_Callback(hObject, eventdata, handles)
% hObject    handle to menu_stackViewer (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


if strcmp(handles.menu_stackViewer.Checked,'on')
    handles.menu_stackViewer.Checked = 'off';
    stackViewer_CloseRequest_Fcn(handles.stackViewer.window, 1, hObject);
else
    handles.menu_stackViewer.Checked = 'on';
    handles.stackViewer.window = figure('Name', 'Stack Viewer', 'MenuBar', 'none', 'Toolbar', 'figure', 'CloseReq', {@stackViewer_CloseRequest_Fcn, hObject}); % 'Keypress', {@onkeypress, hObject}'ResizeFcn', {@movie_slider_ResizeFcn, hObject});
     
    handles.stackViewer.axes = axes('Hittest', 'off');
 
    handles.stackViewer.z_slider = uicontrol(gcf, 'Style', 'slider', 'Position', [115, 25, 350, 20], 'Value', 1, 'Callback', {@stackViewer_z_slider_Callback, hObject});
    uicontrol('Style', 'text', 'String', 'Z', 'Position', [85, 25, 10, 20])
    handles.stackViewer.z_value = 1;
    handles.stackViewer.c_slider = uicontrol(gcf, 'Style', 'slider', 'Position', [115, 5, 350, 20], 'Value', 1, 'Callback', {@stackViewer_c_slider_Callback, hObject});
    uicontrol('Style', 'text', 'String', 'C', 'Position', [85, 5, 10, 20])
    
    uicontrol('Units', 'normalized', 'Style', 'text', 'String', 'Intensity', 'Position', [0.035 0.55 0.1 0.05]); 
    uicontrol('Units', 'normalized', 'Style', 'text', 'String', 'Min:', 'Position', [0.02 0.5 0.05 0.05]); 
    handles.stackViewer.min_box = uicontrol('Units', 'normalized', 'Style', 'edit', 'Position', [0.08 0.505 0.05 0.05], 'Callback', {@stackViewer_minbox_Callback, hObject});
    uicontrol('Units', 'normalized', 'Style', 'text', 'String', 'Max:', 'Position', [0.02 0.4 0.05 0.05]);
    handles.stackViewer.max_box = uicontrol('Units', 'normalized', 'Style', 'edit', 'Position', [0.08 0.405 0.05 0.05], 'Callback', {@stackViewer_maxbox_Callback, hObject});
    
    uicontrol('Units', 'normalized', 'Style', 'text', 'String', 'Show dots', 'Position', [0.01 0.3 0.15 0.05]);
    handles.stackViewer.showdots_box = uicontrol('Units', 'normalized', 'Style', 'checkbox', 'Value', 0, 'Position', [0.14 0.32 0.03 0.03], 'Callback', {@stackViewer_showdots_box_Callback, hObject});   
    
    handles.stackViewer.c_value = 1;
    handles.stackViewer.z_slider.Units = 'normalized';
    handles.stackViewer.c_slider.Units = 'normalized';
    handles = update_stackViewer_display(hObject,handles);   
end
guidata(hObject, handles)

% --------------------------------------------------------------------
function stackViewer_showdots_box_Callback(parent,~ ,hObject)

   handles = guidata(hObject);
   handles = update_stackViewer_display(hObject, handles);
   guidata(hObject, handles);

% --------------------------------------------------------------------
function stackViewer_CloseRequest_Fcn(parent,~ ,hObject)
delete(parent);
handles = guidata(hObject);
handles.menu_stackViewer.Checked = 'off';
guidata(hObject, handles)

% --------------------------------------------------------------------
function stackViewer_z_slider_Callback(parent,~ ,hObject)
handles = guidata(hObject);
z_value = round(get(parent,'Value'));
set(parent,'Value',z_value);

if z_value ~= handles.stackViewer.z_value % New image
    handles.stackViewer.z_value = get(parent,'Value');
    
    handles = update_stackViewer_display(hObject,handles);
    
end


guidata(hObject, handles)

% --------------------------------------------------------------------
function stackViewer_c_slider_Callback(parent,~ ,hObject)
handles = guidata(hObject);

c_value = round(get(parent,'Value'));
set(parent,'Value',c_value);

if c_value ~= handles.stackViewer.c_value % New image
    handles.stackViewer.c_value = get(parent,'Value');
    
    handles = update_stackViewer_display(hObject,handles);
    
end

guidata(hObject, handles)

% --------------------------------------------------------------------
function stackViewer_minbox_Callback(~, ~, hObject)

    handles = guidata(hObject);
    channel = handles.stackViewer.c_value;
    handles.stackViewer.min_value(channel) = str2double(handles.stackViewer.min_box.String);
    handles = update_stackViewer_display(hObject,handles);
    guidata(hObject, handles);

% --------------------------------------------------------------------
function stackViewer_maxbox_Callback(~, ~, hObject)

    handles = guidata(hObject);
    channel = handles.stackViewer.c_value;
    handles.stackViewer.max_value(channel) = str2double(handles.stackViewer.max_box.String);
    handles = update_stackViewer_display(hObject,handles);
    guidata(hObject, handles);

% --------------------------------------------------------------------
function handles = update_display(hObject,handles)
% TO DO: show cell outline
% TO DO: select cell to zoom on, keep focus on
if handles.image_displayed == 0 && handles.stack.number_images >= 1
    handles.image_displayed = 1;
    handles.main_window_axes = gca;
end

if handles.image_displayed
    set(handles.txt_image_counter,'String', [num2str(handles.image_displayed) '/' num2str(handles.stack.number_images) ], 'FontSize', 12);
    if handles.stack.number_images > 1
        set(handles.slider_select_image, 'Enable', 'on');
        set(handles.slider_select_image, 'Value', handles.image_displayed);
        set(handles.slider_select_image, 'Min', 1);
        set(handles.slider_select_image, 'Max', handles.stack.number_images);
        if handles.stack.number_images >1
            set(handles.slider_select_image, 'SliderStep', [1/(handles.stack.number_images-1) , 1/(handles.stack.number_images-1) ]);
        else
            set(handles.slider_select_image, 'SliderStep', [1, 1 ]);
        end
    else
        set(handles.slider_select_image, 'Enable', 'off');
    end
    
    if isempty(handles.current_image_stack)
        % Load the file stack      
        [handles.current_image_stack, handles.num_channels] =  load_position(handles, handles.image_displayed);
    end
    
    num_channels = handles.num_channels;
    imdata = handles.current_image_stack;
    % Show max proj of current stack (first 3 channels)
    max_proj={};max_proj_rescaled=zeros(size(imdata{1},1),size(imdata{1},2),min(num_channels,3));
    for i=1:min(num_channels,3)
        max_proj{i} = max(imdata{i},[],3);
        
        % Automatic rescaling to display RGB
        bincounts = histc(max_proj{i}(:), [0.5:1:65535.5]);
        cdf=cumsum(bincounts)/size(max_proj{i},1)/size(max_proj{i},2);
        
        lower_bound(i) = find(cdf> 0.3,1, 'first');
        higher_bound(i) = find(cdf> 0.999,1, 'first');
        max_proj_rescaled(:,:,i)=mat2gray(max_proj{i},[lower_bound(i) higher_bound(i)]);
        
    end
    
    mask_im = max_proj_rescaled;
    if strcmp(handles.menu_showSegment.Checked, 'on')
        colors = parula(handles.stack.frame(handles.image_displayed).number_cells);
        for k = 1:handles.stack.frame(handles.image_displayed).number_cells
            dilated_im = imdilate(handles.stack.frame(handles.image_displayed).cell(k).mask, strel('square', 5));
            outline_im = imdilate(dilated_im, strel('square', 5)) - handles.stack.frame(handles.image_displayed).cell(k).mask;
            mask_im = imoverlay(mask_im, logical(outline_im), colors(k, :));
        end
    end
        
    axes(handles.main_window_axes);
    imshow(mask_im)
        
    if strcmp(get(handles.menu_stackViewer,'checked'),'on') %also update stackViewer
        handles = update_stackViewer_display(hObject,handles);
    end
    
    guidata(hObject, handles);
     
end

% --------------------------------------------------------------------
function handles = update_stackViewer_display(hObject,handles)

if handles.image_displayed
    
     if isempty(handles.current_image_stack)
        % Load the file stack       
        [handles.current_image_stack, handles.num_channels] = load_position(handles, handles.image_displayed);
        
    end
    
    set(handles.stackViewer.c_slider, 'Min', 1);
    set(handles.stackViewer.c_slider, 'Max', handles.num_channels);
    if handles.num_channels >1
        set(handles.stackViewer.c_slider, 'SliderStep', [1/(handles.num_channels-1) , 1/(handles.num_channels-1) ]);
    else
        set(handles.stackViewer.c_slider, 'SliderStep', [1 , 1]);   
    end
    
    current_channel = handles.stackViewer.c_value;
    handles.num_z_slices = size(handles.current_image_stack{current_channel},3);
    
    set(handles.stackViewer.z_slider, 'Min', 1);
    set(handles.stackViewer.z_slider, 'Max',  handles.num_z_slices );
    
    if handles.num_z_slices >1
        set(handles.stackViewer.z_slider, 'SliderStep', [1/( handles.num_z_slices -1) , 1/( handles.num_z_slices -1) ]);
        set(handles.stackViewer.z_slider, 'Enable', 'on');
    else
        set(handles.stackViewer.z_slider, 'SliderStep', [1 , 1]);
        set(handles.stackViewer.z_slider, 'Enable', 'off');
        
    end
        
    current_z_slice = handles.stackViewer.z_value;
    
    if current_z_slice > handles.num_z_slices
        current_z_slice = handles.num_z_slices;
        handles.stackViewer.z_value = current_z_slice;
        set(handles.stackViewer.z_slider, 'Value', current_z_slice);
    end
    
    imdata = handles.current_image_stack;

    current_image = imdata{current_channel}(:,:,current_z_slice);
   
    % Automatic rescaling
    bincounts = histc(current_image(:), [0.5:1:65535.5]);
    cdf=cumsum(bincounts)/size(current_image,1)/size(current_image,2);
    
    if ~isfield(handles.stackViewer, 'min_value')
        handles.stackViewer.min_value = zeros(1, handles.num_channels);
        handles.stackViewer.max_value = zeros(1, handles.num_channels);
    end
    
    if handles.stackViewer.min_value(current_channel) == 0
        lower_bound = find(cdf> 0.01,1, 'first');
        handles.stackViewer.min_value(current_channel) = lower_bound;
        higher_bound = find(cdf> 0.99,1, 'first');
        handles.stackViewer.max_value(current_channel) = higher_bound;
    end
    
    set(handles.stackViewer.min_box, 'String', num2str(handles.stackViewer.min_value(current_channel)));
    set(handles.stackViewer.max_box, 'String', num2str(handles.stackViewer.max_value(current_channel)));
    
    figure(handles.stackViewer.window);
    xlim = get(gca, 'XLim');
    ylim = get(gca, 'YLim');
    
    mask_im = mat2gray(current_image,[handles.stackViewer.min_value(current_channel) handles.stackViewer.max_value(current_channel)]);
    if strcmp(handles.menu_showSegment.Checked, 'on')
        colors = parula(handles.stack.frame(handles.image_displayed).number_cells);
        for k = 1:handles.stack.frame(handles.image_displayed).number_cells
            dilated_im = imdilate(handles.stack.frame(handles.image_displayed).cell(k).mask, strel('square', 5));
            outline_im = imdilate(dilated_im, strel('square', 5)) - handles.stack.frame(handles.image_displayed).cell(k).mask;
            mask_im = imoverlay(mask_im, logical(outline_im), colors(k, :));
        end
    end
    
    imshow(mask_im)
    zoom(gcf,'reset');
    if sum(xlim ~= [0 1]) && sum(ylim ~= [0 1]) % Already something displayed
        set(gca,'xlim',xlim);
        set(gca,'ylim',ylim);
    end
    
    if isfield(handles.stack.frame(handles.image_displayed), 'cell') && ~isempty(handles.stack.frame(handles.image_displayed).cell) 
        if get(handles.stackViewer.showdots_box, 'Value') && isfield(handles.stack.frame(handles.image_displayed).cell(1), 'dots')
            for cellNo=1:length(handles.stack.frame(handles.image_displayed).cell)
                if ~isempty(handles.stack.frame(handles.image_displayed).cell(cellNo).dots)
                    good_dots_stats = handles.stack.frame(handles.image_displayed).cell(cellNo).dots(current_channel).properties;
                    centroid_to_plot = zeros(length(good_dots_stats),3);
                    
                    if handles.num_z_slices ==1
                        for p = 1:length(good_dots_stats)
                            centroid_to_plot(p,:) = [good_dots_stats(p).Centroid 1];
                        end
                    else
                        for p = 1:length(good_dots_stats)
                            centroid_to_plot(p,:) = good_dots_stats(p).Centroid;
                        end
                    end


                    dots_in_plane = round(centroid_to_plot(:,3)) == current_z_slice;
                    centroids_in_plane = centroid_to_plot(dots_in_plane,1:2);
                    hold on;
                    viscircles(centroid_to_plot(:,1:2),repmat(2,size(centroid_to_plot(:,1:2),1),1),'Edgecolor','g');
                    viscircles(centroids_in_plane,repmat(2,size(centroids_in_plane,1),1));
                    hold off;
                end
                
            end 
        end
    end
    
end

% --------------------------------------------------------------------
function menu_addCell_Callback(hObject, eventdata, handles)
% hObject    handle to menu_addCell (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if handles.image_displayed
    h = figure;
    imdata = handles.current_image_stack;
    % Show max proj of current stack (first 3 channels)
    max_proj={};max_proj_rescaled=zeros(size(imdata{1},1),size(imdata{1},2),min(handles.num_channels,3));
    for i=1:min(handles.num_channels,3)
        max_proj{i} = max(imdata{i},[],3);
        
        % Automatic rescaling to display RGB
        bincounts = histc(max_proj{i}(:), [0.5:1:65535.5]);
        cdf=cumsum(bincounts)/size(max_proj{i},1)/size(max_proj{i},2);
        
        lower_bound(i) = find(cdf> 0.3,1, 'first');
        higher_bound(i) = find(cdf> 0.999,1, 'first');
        max_proj_rescaled(:,:,i)=mat2gray(max_proj{i},[lower_bound(i) higher_bound(i)]);
        
    end
    %          figure(handles.ax_main_image);
    imshow(max_proj_rescaled)
    
    BW_mask = roipoly;
    handles.stack.frame(handles.image_displayed).number_cells = handles.stack.frame(handles.image_displayed).number_cells + 1;

    
    if isfield(handles.stack, 'dot_thresholds') && handles.stack.live_dot_finding
        [handles, dots] = dotCounting(handles, imdata, BW_mask);
        handles.stack.frame(handles.image_displayed).cell(handles.stack.frame(handles.image_displayed).number_cells).dots = dots;
    end
    
    close(h);
    set(handles.menu_removeCell, 'Enable', 'On')
    
    handles.stack.frame(handles.image_displayed).cell(handles.stack.frame(handles.image_displayed).number_cells).mask  = BW_mask;
    
    handles = update_display(hObject,handles);
    guidata(hObject, handles)  
end

% --------------------------------------------------------------------
function menu_segmentCells_Callback(hObject, eventdata, handles)
% hObject    handle to menu_segmentCells (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    if strcmp(hObject.Checked, 'off')
        set(hObject, 'Checked', 'on')
    else
        figure(handles.seg_parameters.window);
    end
    
    handles.seg_parameters.window = figure('Name', 'Segmentation Parameters', 'Menubar', 'none', 'Toolbar', 'figure', 'Position', [500 500 200 100], 'CloseRequestFcn', {@menu_segmentCells_CloseRequestFcn, hObject});
    ax = axes('Units', 'pixels', 'Visible', 'off'); 
    uicontrol('Style', 'text', 'String', 'Channel', 'Position', [5 60 50 10]);  
    uicontrol('Style', 'text', 'String', 'Threshold', 'Position', [5 40 50 10]); 
    uicontrol('Style', 'text', 'String', 'Channel', 'Position', [115 60 50 10]); 
    uicontrol('Style', 'text', 'String', 'Threshold', 'Position', [115 40 50 10]);
    uicontrol('Style', 'text', 'String', 'Cell', 'Position', [30 80 50 10]);
    uicontrol('Style', 'text', 'String', 'Nucleus', 'Position', [135 80 50 10]);
    
    if ~isfield(handles.seg_parameters, 'seg_channel')
        handles.seg_parameters.seg_channel = 1; handles.seg_parameters.seg_threshold = 100; 
        handles.seg_parameters.nuc_channel = 1; handles.seg_parameters.nuc_threshold = 1000;
    end
    
    handles.seg_parameters.seg_channel_editbox = uicontrol('Style', 'edit', 'String', num2str(handles.seg_parameters.seg_channel), 'Position', [55 55 30 20], 'Callback', {@seg_channel_Callback, hObject});
    handles.seg_parameters.seg_threshold_editbox = uicontrol('Style', 'edit', 'String', num2str(handles.seg_parameters.seg_threshold), 'Position', [55 35 30 20], 'Callback', {@seg_threshold_Callback, hObject});
    handles.seg_parameters.nuc_channel_editbox = uicontrol('Style', 'edit', 'String', num2str(handles.seg_parameters.nuc_channel), 'Position', [165 55 30 20], 'Callback', {@nuc_channel_Callback, hObject});
    handles.seg_parameters.nuc_threshold_editbox = uicontrol('Style', 'edit', 'String', num2str(handles.seg_parameters.nuc_threshold), 'Position', [165 35 30 20], 'Callback', {@nuc_threshold_Callback, hObject});
    
    uicontrol('Style', 'pushbutton', 'String', 'Segment!', 'Position', [80 5 55 20], 'Callback', {@seg_button_Callback, hObject});
    guidata(hObject, handles);

% --------------------------------------------------------------------
function menu_segmentCells_CloseRequestFcn(parent, ~, hObject)
    
    handles = guidata(hObject);
    set(handles.menu_segmentCells, 'Checked', 'off')
    delete(parent);
    guidata(hObject, handles);
    
% --------------------------------------------------------------------   
function seg_channel_Callback(parent, ~, hObject)
        
    handles = guidata(hObject);
    handles.seg_parameters.seg_channel = str2double(get(handles.seg_parameters.seg_channel_editbox, 'String'));
    guidata(hObject, handles);
    
% --------------------------------------------------------------------
function seg_threshold_Callback(~, ~, hObject)

    handles = guidata(hObject);
    handles.seg_parameters.seg_threshold = str2double(get(handles.seg_parameters.seg_threshold_editbox, 'String'));
    
    imdata =  handles.current_image_stack;
    seg_channel_data = imdata{handles.seg_parameters.seg_channel};
    nuc_channel_data = imdata{handles.seg_parameters.nuc_channel};
    [seg_im, ~] = segment_on_bg(seg_channel_data, handles.seg_parameters.seg_threshold, nuc_channel_data, handles.seg_parameters.nuc_threshold);
    
    for p = 1:max(max(seg_im))
        BW_mask = 0*seg_im; BW_mask(seg_im == p) = 1;
        handles.stack.frame(handles.image_displayed).number_cells = handles.stack.frame(handles.image_displayed).number_cells + 1;
        handles.stack.frame(handles.image_displayed).cell(handles.stack.frame(handles.image_displayed).number_cells).mask = BW_mask;
    end
    handles = update_display(hObject,handles);
    guidata(hObject, handles);

% --------------------------------------------------------------------
function nuc_channel_Callback(~, ~, hObject)
        
    handles = guidata(hObject);
    handles.seg_parameters.nuc_channel = str2double(get(handles.seg_parameters.nuc_channel_editbox, 'String'));
    guidata(hObject, handles);
    
% --------------------------------------------------------------------
function nuc_threshold_Callback(~, ~, hObject)

    handles = guidata(hObject);
    handles.seg_parameters.nuc_threshold = str2double(get(handles.seg_parameters.nuc_threshold_editbox, 'String'));
    imdata =  handles.current_image_stack;
    seg_channel_data = imdata{handles.seg_parameters.seg_channel};
    nuc_channel_data = imdata{handles.seg_parameters.nuc_channel};
    [seg_im, ~] = segment_on_bg(seg_channel_data, handles.seg_parameters.seg_threshold, nuc_channel_data, handles.seg_parameters.nuc_threshold);
    
    for p = 1:max(max(seg_im))
        BW_mask = 0*seg_im; BW_mask(seg_im == p) = 1;
        handles.stack.frame(handles.image_displayed).number_cells = handles.stack.frame(handles.image_displayed).number_cells + 1;
        handles.stack.frame(handles.image_displayed).cell(handles.stack.frame(handles.image_displayed).number_cells).mask = BW_mask;
    end
    handles = update_display(hObject,handles);
    guidata(hObject, handles);

% --------------------------------------------------------------------
function seg_button_Callback(~, ~, hObject)
     
    handles = guidata(hObject);
    % load each image and segment
    for k = 1:handles.stack.number_images
        imdata = load_position(handles,k);
        seg_channel_data = imdata{handles.seg_parameters.seg_channel};
        nuc_channel_data = imdata{handles.seg_parameters.nuc_channel};
        [seg_im, ~] = segment_on_bg(seg_channel_data, handles.seg_parameters.seg_threshold, nuc_channel_data, handles.seg_parameters.nuc_threshold);
        
        for p = 1:max(max(seg_im))
            BW_mask = 0*seg_im; BW_mask(seg_im == p) = 1;
            handles.stack.frame(k).number_cells = handles.stack.frame(k).number_cells + 1;
            handles.stack.frame(k).cell(handles.stack.frame(k).number_cells).mask = BW_mask;
        end
    end
    handles = update_display(hObject,handles);
    guidata(hObject, handles);

% --------------------------------------------------------------------
function menu_showSegment_Callback(hObject, eventdata, handles)
% hObject    handle to menu_showSegment (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    
    if strcmp(get(hObject, 'Checked'), 'off')
        set(hObject, 'Checked', 'on');
    else
        set(hObject, 'Checked', 'off');
    end    
    handles = update_display(hObject,handles);
    guidata(hObject, handles);

% --------------------------------------------------------------------
function menu_detectDots_Callback(hObject, eventdata, handles)
% hObject    handle to menu_detectDots (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    if strcmp(hObject.Checked,'on')
%         set(hObject, 'Checked', 'off')
        figure(handles.threshold.window);
    else
        hObject.Checked = 'on';
        window_position = getpixelposition(handles.figure1);
        handles.threshold.window = figure('Name', 'Threshold settings', 'Menubar', 'none', 'Position', [window_position(1), window_position(2), 261, 271], 'CloseReq', {@menu_threshold_CloseRequest_Fcn, hObject}); % 'Keypress', {@onkeypress, hObject}'ResizeFcn', {@movie_slider_ResizeFcn, hObject});
        handles.threshold.status = uicontrol(gcf, 'Style', 'text', 'Position', [141, 150, 104, 35], 'String', 'Status: Done', 'FontSize', 16);
        
        if ~isfield(handles.stack, 'dot_thresholds')
            handles = get_thresholds(handles);
            if handles.stack.frame(handles.image_displayed).number_cells > 0
                handles = rethreshold_cell(handles, handles.image_displayed, length(handles.stack.frame(handles.image_displayed).cell));
            end
        end
        
        handles.threshold.chk_auto_threshold = uicontrol(gcf, 'Style', 'checkbox', 'Position', [31, 218, 95, 19], 'Value', 1, 'String', 'Auto Threshold', 'Callback', {@chk_auto_threshold_Callback, hObject});
       
        if isfield(handles.stack,'autothreshold')
            set(handles.threshold.chk_auto_threshold, 'value',handles.stack.autothreshold);
        else
            handles.stack.autothreshold = 1;
        end
        handles.threshold.chk_live_dot_finding = uicontrol(gcf, 'Style', 'checkbox', 'Position', [31, 250, 95, 19], 'Value', 1, 'String', 'Live dot finding', 'Callback', {@chk_live_dot_finding_Callback, hObject});
        
        if isfield(handles.stack,'live_dot_finding')
            set(handles.threshold.chk_live_dot_finding, 'value',handles.stack.live_dot_finding);
        else
            handles.stack.live_dot_finding = 1;
        end
      
    
        handles.threshold.button_rethreshold = uicontrol(gcf, 'Style', 'pushbutton', 'Position', [141, 210, 84, 35], 'String', 'Detect Dots', 'Callback', {@button_rethreshold_Callback, hObject});
        
        for i=1:handles.num_channels
            handles.threshold.(['edit_channel' num2str(i)]) = uicontrol(gcf, 'Style', 'edit', 'Position', [26, 170-40*(i-1), 46, 34], 'String', num2str(handles.stack.dot_thresholds(i)), 'Callback', {@change_threshold_Callback, hObject}, 'Tag', ['edit_channel' num2str(i)], 'Enable', 'on');
            %handles.threshold.(['channel' num2str(i)]) = '';
            if handles.stack.autothreshold
                set(handles.threshold.(['edit_channel' num2str(i)]),'enable','off');
            end
        end
        
    end
    guidata(hObject, handles)
    
% --------------------------------------------------------------------
   
function chk_live_dot_finding_Callback(parent,~ ,hObject)

handles = guidata(hObject);

if get(parent,'value')
    handles.stack.live_dot_finding = 1;
else
    handles.stack.live_dot_finding = 0;
end

guidata(hObject, handles)
% --------------------------------------------------------------------
function handles = get_thresholds(handles)
    
    [imdata, num_channels] = load_position(handles,1);

    seg_im = 0*imdata{1}(:, :, 1);
    if handles.stack.frame(1).number_cells > 0
        for p = 1:handles.stack.frame(1).number_cells
            seg_im(handles.stack.frame(1).cell(p).mask == 1) = p;
        end
    end

    handles.stack.dot_thresholds = detect_dot_thresh(imdata, num_channels, seg_im);
    
      
% --------------------------------------------------------------------
function menu_removeCell_Callback(hObject, eventdata, handles)
% hObject    handle to menu_removeCell (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    if ~isempty(handles.stack.frame(handles.image_displayed).cell)
        h= figure;
        num_channels = handles.num_channels;
        imdata = handles.current_image_stack;
        % Show max proj of current stack (first 3 channels)
        max_proj={};max_proj_rescaled=zeros(size(imdata{1},1),size(imdata{1},2),min(num_channels,3));
        for i=1:min(num_channels,3)
            max_proj{i} = max(imdata{i},[],3);

            % Automatic rescaling to display RGB
            bincounts = histc(max_proj{i}(:), [0.5:1:65535.5]);
            cdf=cumsum(bincounts)/size(max_proj{i},1)/size(max_proj{i},2);

            lower_bound(i) = find(cdf> 0.3,1, 'first');
            higher_bound(i) = find(cdf> 0.999,1, 'first');
            max_proj_rescaled(:,:,i)=mat2gray(max_proj{i},[lower_bound(i) higher_bound(i)]);

        end
        %          figure(handles.ax_main_image);

        mask_im = max_proj_rescaled;
        colors = parula(handles.stack.frame(handles.image_displayed).number_cells);
        for k = 1:handles.stack.frame(handles.image_displayed).number_cells
            dilated_im = imdilate(handles.stack.frame(handles.image_displayed).cell(k).mask, strel('square', 5));
            outline_im = dilated_im - handles.stack.frame(handles.image_displayed).cell(k).mask;
            mask_im = imoverlay(mask_im, logical(outline_im), colors(k, :));
        end

        imshow(mask_im)


        [xClick, yClick, W]=ginput(1);
        xClick = round(xClick);
        yClick = round(yClick);
        click_image = zeros(size(max_proj_rescaled)) ;
        click_image(yClick,xClick) = 1;
        cellToRemove = 0;
        for i=1:handles.stack.frame(handles.image_displayed).number_cells
            match = find(logical(click_image(:, :, 1)).*handles.stack.frame(handles.image_displayed).cell(i).mask,1,'first');
            if match
                cellToRemove = i;
            end
        end
        if cellToRemove
            indexToKeep = setxor(1:handles.stack.frame(handles.image_displayed).number_cells,cellToRemove);
            handles.stack.frame(handles.image_displayed).number_cells = handles.stack.frame(handles.image_displayed).number_cells -1;
            if handles.stack.frame(handles.image_displayed).number_cells == 0
                set(handles.menu_removeCell, 'Enable', 'On');
            end

            handles.stack.frame(handles.image_displayed).cell = handles.stack.frame(handles.image_displayed).cell(indexToKeep);
        end
        close(h);
    end
    handles = update_display(hObject,handles);
    guidata(hObject, handles)

% --------------------------------------------------------------------
function menu_newStack_Callback(hObject, eventdata, handles)
% hObject    handle to menu_newStack (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


if handles.open_stack == 0
    [FileName,PathName,FilterIndex] = uiputfile('*.mat');
    if FileName
        handles.stack.name = [PathName FileName];
        handles.open_stack = 1;
        set(handles.txt_stack,'String', FileName)
        handles.stack.number_images = 0;
        handles.stack.image_path_cell = {};
        handles.stack.frame.cell = struct;
        handles.stack.image_path_number = [];
        handles.stack.frame.number_cells = 0;
        handles.stack.live_dot_finding = 1;
        handles.image_displayed = 0;
        handles.current_image_stack = [];
        
    end
else
    msgbox('A stack is already open. Close current stack')
    
end


guidata(hObject, handles)

% TO DO : save threshold settings in stack and load it
% TO DO : close stack

% --------------------------------------------------------------------
function menu_openStack_Callback(hObject, eventdata, handles)
% hObject    handle to menu_openStack (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if handles.open_stack == 0
    [FileName,PathName] = uigetfile('*.mat');
    if FileName
        load([PathName FileName]);
        handles.stack = stack;
        handles.open_stack = 1;
        set(handles.txt_stack,'String', FileName)
        handles.image_displayed = 0;
        handles.current_image_stack = [];
        handles = update_display(hObject,handles);
        
    end
else
    msgbox('A stack is already open. Close current stack')
    
end
    guidata(hObject, handles)

% --------------------------------------------------------------------
function menu_saveStack_Callback(hObject, eventdata, handles)
% hObject    handle to menu_saveStack (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    stack = handles.stack;
    save(handles.stack.name, 'stack')


% --------------------------------------------------------------------
function menu_saveStackAs_Callback(hObject, eventdata, handles)
% hObject    handle to menu_saveStackAs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    
    stack = handles.stack;
    [FileName,PathName] = uiputfile('*.mat');
    if FileName
        set(handles.txt_stack,'String', FileName)
        handles.stack.name = [PathName FileName];
        save(handles.stack.name, 'stack');
    end
    guidata(hObject, handles);
    
    
% -------------------CZI images----------------------------------------
function menu_addImages_Callback(hObject, eventdata, handles)
% hObject    handle to menu_addImages (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    [FileName,PathName] = uigetfile('*.czi','MultiSelect','on');

    if ~isfield(handles.stack, 'name')
            [StackFileName,StackPathName] = uiputfile('*.mat', 'Specify Stack name');
            handles.stack.name = [StackPathName, StackFileName];
            handles.open_stack = 1;
            set(handles.txt_stack,'String', StackFileName)
            
            handles.stack.number_images = 0;
            handles.stack.image_path_cell = {};
            handles.stack.frame.cells = struct;
            handles.stack.frame.number_cells = 0;
            handles.image_displayed = 0;
            handles.current_image_stack = [];
            
            guidata(hObject, handles)
    end


    if ~isa(FileName, 'double') %Dialog not closed

        if isa(FileName, 'cell') % More than one file

            for i=1:length(FileName)
                handles.stack.image_path_cell{length( handles.stack.image_path_cell) + 1} = [PathName FileName{i}];
                handles.stack.frame(handles.stack.number_images + i).number_cells = 0;
            end

            handles.stack.number_images = handles.stack.number_images + length(FileName);


        else
            if isfield(handles.stack, 'number_images')
                handles.stack.number_images = handles.stack.number_images + 1;
            else
                handles.stack.number_images = 1;
            end

            handles.stack.frame(handles.stack.number_images).number_cells = 0;
            handles.stack.image_path_cell{length( handles.stack.image_path_cell) + 1} = [PathName FileName];
        end
        handles.stack.filetype = 'czi';
        handles = update_display(hObject,handles);
        guidata(hObject, handles)


    end
    % --------------------------------------------------------------------
function menu_addImagesND_Callback(hObject, eventdata, handles)
% hObject    handle to menu_addImagesND (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[FileName,PathName] = uigetfile('*.tif','Pick one file with basename');
[startIndex,endIndex] = regexp(FileName,'_w.*_s');
basename =FileName(1:(startIndex-1));
if ~isa(FileName, 'double') %Dialog not closed
    
    files = dir([PathName filesep basename '*.tif']);
    number_position = 1;
    number_channel = 1;
    for i=1:length(files)
        FileName=files(i).name;
        [startIndex,endIndex] = regexp(FileName,'_w.*_s');
        channel_num = str2num(FileName(startIndex+2));
        [startIndex,endIndex] = regexp(FileName,'_s.*(?i:tif)');
        pos_num = str2num(FileName(startIndex+2:endIndex-4));
        if pos_num > number_position
            number_position = pos_num;
        end
        if channel_num > number_channel
            number_channel = channel_num;
        end
    end
    prompt = ['Detected ' num2str(number_position) ' positions with ' num2str(number_channel) ' channels\n' 'Add positions (e.g. 1-5):'];
    dlg_title = 'Select positions';
    answer = inputdlg(prompt,dlg_title);
    if ~isempty(answer)
        %% PUT BACK
        index = regexp(answer{1}, '-');
        pos_start = str2num(answer{1}(1:index - 1));
        pos_end = str2num(answer{1}(index + 1:end));
        if logical(pos_start) && logical(pos_end)
            pos_vector = pos_start:pos_end;
            for i=1:length(pos_vector)
                handles.stack.image_path_cell{length( handles.stack.image_path_cell) + 1} = [PathName basename];
                handles.stack.image_path_number(length( handles.stack.image_path_number) + 1) = pos_vector(i);
                handles.stack.frame(handles.stack.number_images + i).number_cells = 0;
            end
            
            handles.stack.number_images = handles.stack.number_images + pos_end-pos_start + 1;
            handles.stack.filetype = 'tif';
            handles = update_display(hObject,handles);
            guidata(hObject, handles)
            
        end
    end
  
end



% --------------------------------------------------------------------
function [imdata, num_channels] = load_position(handles, position)

switch handles.stack.filetype
    case 'czi'
        [imdata, num_channels] = czi_open(handles.stack.image_path_cell{position});
    case 'tif'        
        [imdata, num_channels] = tif_open(handles.stack.image_path_cell{position},handles.stack.image_path_number(position));       
end
% --------------------------------------------------------------------

function [imdata, num_channels] = tif_open(basename_path,actual_position)
files = dir([basename_path '*s' num2str(actual_position) '.TIF']);
num_channels=length(files);
imdata={};
for i=1:num_channels
    %% PUT BACK
    imdata{i} = loadStack([files(i).name]);
end


    % --------------------------------------------------------------------
function menu_removeImage_Callback(hObject, eventdata, handles)
% hObject    handle to menu_removeImage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


    if handles.stack.number_images > 1
        handles.stack.number_images = handles.stack.number_images - 1;
        handles.stack.image_path_cell = handles.stack.image_path_cell([1:handles.image_displayed - 1, handles.image_displayed + 1:end]);
        handles.stack.frame = handles.stack.frame([1:handles.image_displayed - 1, handles.image_displayed + 1:end]);
        handles.image_displayed = max(1, handles.image_displayed - 1);
        handles.current_image_stack = [];
        handles = update_display(hObject,handles);
        guidata(hObject, handles)
    end

% --------------------------------------------------------------------
function [handles, dots] = dotCounting(handles, imdata, cell_mask)
    set(gcf,'pointer','watch');

    BW_mask = cell_mask;
    dots = detect_dots(imdata, BW_mask, handles.num_channels, handles.stack.dot_thresholds);
    
    % TO DO
    % -remove dots that colocalize
    set(gcf,'pointer','arrow');

% --- Executes on slider movement.
function slider_select_image_Callback(hObject, eventdata, handles)
% hObject    handle to slider_select_image (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

    image_to_display = round(get(hObject,'Value'));
    set(hObject,'Value',image_to_display);

    if image_to_display ~= handles.image_displayed % New image
        handles.image_displayed = get(hObject,'Value');
        handles.current_image_stack = [];


        handles = update_display(hObject,handles);

    end

    guidata(hObject, handles)


% --- Executes during object creation, after setting all properties.
function slider_select_image_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_select_image (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --------------------------------------------------------------------
function chk_auto_threshold_Callback(parent,~ ,hObject)

    handles = guidata(hObject);

    if get(parent,'value')
        handles.stack.autothreshold = 1;
        handles.threshold.status.String = 'Detecting ...';
         pause(1);
        handles = get_thresholds(handles);
        if handles.stack.frame(handles.image_displayed).number_cells > 0
            handles = rethreshold_cell(handles, handles.image_displayed, length(handles.stack.frame(handles.image_displayed).cell));
            if strcmp(get(handles.menu_stackViewer,'checked'),'on') %also update display
                handles = update_stackViewer_display(hObject,handles);
            end
        end
        
        for i=1:handles.num_channels
            handles.threshold.(['edit_channel' num2str(i)]).String = num2str(handles.stack.dot_thresholds(i));
            handles.threshold.(['edit_channel' num2str(i)]).Enable = 'off';
        end
        handles.threshold.status.String = 'Status: Done';
    else
        handles.stack.autothreshold = 0;
        for i=1:handles.num_channels
            handles.threshold.(['edit_channel' num2str(i)]).Enable = 'on';
        end
    end

    guidata(hObject, handles)

% --------------------------------------------------------------------
function button_rethreshold_Callback(parent,~ ,hObject)

    handles = guidata(hObject);
    handles.threshold.status.String = 'Detecting ...';
    pause(1);
    h = waitbar(0);
    mean_detection_time = 0;
    tic
    for fov=1:handles.stack.number_images
       if  ~isempty(handles.stack.frame(fov).cell)
           for cellNo=1:handles.stack.frame(fov).number_cells
               handles = rethreshold_cell(handles, fov, cellNo);
               mean_detection_time = toc/fov;
               % TO DO: threshold all cells in fov at same time
           end
       end
       waitbar(fov/handles.stack.number_images, h, ['Time remaining: ', num2str(mean_detection_time*(handles.stack.number_images - fov)), ' secs']);
    end
    close(h);
    handles.threshold.status.String = 'Status: Done';
    guidata(hObject, handles)

% --------------------------------------------------------------------
function change_threshold_Callback(parent,~ ,hObject)
    handles = guidata(hObject);
    parent_tag = get(parent,'Tag');
    channel_to_change = str2double(parent_tag(end));
    old_threshold = handles.stack.dot_thresholds(channel_to_change);
    new_threshold = str2double(get(parent,'String'));

    if ~handles.stack.autothreshold && ~isempty(new_threshold) %input validation
        handles.stack.dot_thresholds(channel_to_change) = new_threshold;
        %handles.threshold.(['channel' channel_to_change]) = new_threshold;
    else %put it back the way it was
        set(parent,'String', num2str(old_threshold));
    end
    
    % Rethreshold last cell
    if handles.stack.frame(handles.image_displayed).number_cells > 0  && handles.stack.live_dot_finding
        handles.threshold.status.String = 'Detecting ...';
        pause(1);    
        handles = rethreshold_cell(handles, handles.image_displayed, length(handles.stack.frame(handles.image_displayed).cell));
        if strcmp(get(handles.menu_stackViewer,'checked'),'on') %also update display
            handles = update_stackViewer_display(hObject,handles);
        end
    end
    handles.threshold.status.String = 'Status: Done';
    guidata(hObject, handles)


function menu_threshold_CloseRequest_Fcn(parent,~ ,hObject)
    delete(parent);
    handles = guidata(hObject);
    handles.menu_detectDots.Checked = 'off';
    guidata(hObject, handles)


function handles = rethreshold_cell(handles, frame_num, cell_num)

    BW_mask = handles.stack.frame(frame_num).cell(cell_num).mask ;
    if frame_num ~= handles.image_displayed %need to load the file
        imdata = load_position(handles,frame_num);
    else
        imdata = handles.current_image_stack;
    end
    [handles,dots] = dotCounting(handles, imdata, BW_mask);
    handles.stack.frame(frame_num).cell(cell_num).dots = dots;


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure

 selection = questdlg('Close dotCounting_gui?',...
                     'Close Request',...
                     'Yes','No','No');
 switch selection
   case 'Yes'
     if strcmp(get(handles.menu_stackViewer,'Checked'),'on')
         stackViewer_CloseRequest_Fcn(handles.stackViewer.window, 1, hObject);
     end
     
     if strcmp(get(handles.menu_detectDots,'Checked'),'on')
         menu_threshold_CloseRequest_Fcn(handles.threshold.window, 1, hObject);
     end
      
     if strcmp(get(handles.menu_segmentCells,'Checked'),'on')
         menu_segmentCells_CloseRequestFcn(handles.seg_parameters.window, 1, hObject);
     end
     delete(hObject);
     
   case 'No'
     return
 end


