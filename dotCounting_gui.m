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

% Last Modified by GUIDE v2.5 12-Dec-2017 18:26:51

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
% handles.num_channel = 3;
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
function Untitled_1_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_1 (see GCBO)
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
    handles.stackViewer.window = figure('Name', 'Stack Viewer', 'Menubar', 'none', 'Toolbar', 'figure', 'CloseReq', {@stackViewer_CloseRequest_Fcn, hObject}); % 'Keypress', {@onkeypress, hObject}'ResizeFcn', {@movie_slider_ResizeFcn, hObject});
    handles.stackViewer.axes = axes('Hittest', 'off');
    
    handles.stackViewer.z_slider = uicontrol(gcf, 'Style', 'slider', 'Position', [65, 5, 450, 30], 'Value', 1, 'Callback', {@stackViewer_z_slider_Callback, hObject});
    handles.stackViewer.z_value = 1;
    handles.stackViewer.c_slider = uicontrol(gcf, 'Style', 'slider', 'Position', [65, 35, 450, 30], 'Value', 1, 'Callback', {@stackViewer_c_slider_Callback, hObject});
    handles.stackViewer.c_value = 1;
    handles.stackViewer.z_slider.Units = 'normalized';
    handles.stackViewer.c_slider.Units = 'normalized';
    handles = update_stackViewer_display(hObject,handles);

    
end
guidata(hObject, handles)

function stackViewer_CloseRequest_Fcn(parent,~ ,hObject)
delete(parent);
handles = guidata(hObject);
handles.menu_stackViewer.Checked = 'off';
guidata(hObject, handles)

function stackViewer_z_slider_Callback(parent,~ ,hObject)
handles = guidata(hObject);
z_value = round(get(parent,'Value'));
set(parent,'Value',z_value);

if z_value ~= handles.stackViewer.z_value % New image
    handles.stackViewer.z_value = get(parent,'Value');
    
    handles = update_stackViewer_display(hObject,handles);
    
end


guidata(hObject, handles)

function stackViewer_c_slider_Callback(parent,~ ,hObject)
handles = guidata(hObject);

c_value = round(get(parent,'Value'));
set(parent,'Value',c_value);

if c_value ~= handles.stackViewer.c_value % New image
    handles.stackViewer.c_value = get(parent,'Value');
    
    handles = update_stackViewer_display(hObject,handles);
    
end

guidata(hObject, handles)

%%%%%%%%%%%%%%%%%%%%%%

function handles = update_stackViewer_display(hObject,handles)

if handles.image_displayed
    
     if isempty(handles.current_image_stack)
        % Load the file stack 
        
        [handles.current_image_stack, handles.num_channel] = load_position(handles, handles.image_displayed);
        
    end
    
    
    %     set(handles.txt_image_counter,'String', [num2str(handles.image_displayed) '/' num2str(handles.stack.number_images) ] );
    %     set(handles.slider_select_image, 'Value', handles.image_displayed);
    set(handles.stackViewer.c_slider, 'Min', 1);
    set(handles.stackViewer.c_slider, 'Max', handles.num_channel);
    if handles.num_channel >1
        set(handles.stackViewer.c_slider, 'SliderStep', [1/(handles.num_channel-1) , 1/(handles.num_channel-1) ]);
    else
        set(handles.stackViewer.c_slider, 'SliderStep', [1 , 1]);   
    end
    
    current_channel = handles.stackViewer.c_value;
   
    handles.num_z_slices = size(handles.current_image_stack{current_channel},3);
    current_z_slice = handles.stackViewer.z_value;
    
    if current_z_slice > handles.num_z_slices
        current_z_slice = handles.num_z_slices;
        handles.stackViewer.z_value = current_z_slice;
        set(handles.stackViewer.z_slider, 'Value', current_z_slice);
    end
    
    set(handles.stackViewer.z_slider, 'Min', 1);
    set(handles.stackViewer.z_slider, 'Max',  handles.num_z_slices );
    if handles.num_z_slices >1
        set(handles.stackViewer.z_slider, 'SliderStep', [1/( handles.num_z_slices -1) , 1/( handles.num_z_slices -1) ]);
        set(handles.stackViewer.z_slider, 'Enable', 'on');
    else
        set(handles.stackViewer.z_slider, 'SliderStep', [1 , 1]);
        set(handles.stackViewer.z_slider, 'Enable', 'off');
        
    end
    

    
    imdata = handles.current_image_stack;

    current_image = imdata{current_channel}(:,:,current_z_slice);
    
    
    % Automatic rescaling
    bincounts = histc(current_image(:), [0.5:1:65535.5]);
    cdf=cumsum(bincounts)/size(current_image,1)/size(current_image,2);
    
    lower_bound = find(cdf> 0.01,1, 'first');
    higher_bound = find(cdf> 0.99,1, 'first');
    figure(handles.stackViewer.window);
%     ax = axes('Parent',fig);
    xlim = get(gca, 'XLim');
    ylim = get(gca, 'YLim');
    imshow(current_image,[lower_bound higher_bound])
    zoom(gcf,'reset');
    if sum(xlim ~= [0 1]) && sum(ylim ~= [0 1]) % Already something displayed
        set(gca,'xlim',xlim);
        set(gca,'ylim',ylim);
    end
    if isfield(handles.stack.frame(handles.image_displayed).cell,'dots') 
        for cellNo=1:length(handles.stack.frame(handles.image_displayed).cell)
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
                   viscircles(centroid_to_plot(:,1:2),repmat(2,size(centroid_to_plot(:,1:2),1),1),'color','g');
                   viscircles(centroids_in_plane,repmat(2,size(centroids_in_plane,1),1));
                   hold off;
                   
                
        end
    end
    
    %          figure(handles.ax_main_image);
end
% guidata(hObject, handles)
% --------------------------------------------------------------------
function menu_addCell_Callback(hObject, eventdata, handles)
% hObject    handle to menu_addCell (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if handles.image_displayed
    h = figure;
    imdata = handles.current_image_stack;
    % Show max proj of current stack (first 3 channels)
    max_proj={};max_proj_rescaled=zeros(size(imdata{1},1),size(imdata{1},2),min(handles.num_channel,3));
    for i=1:min(handles.num_channel,3)
        max_proj{i} = max(imdata{i},[],3);
        
        % Automatic rescaling to display RGB
        bincounts = histc(max_proj{i}(:), [0.5:1:65535.5]);
        cdf=cumsum(bincounts)/size(max_proj{i},1)/size(max_proj{i},2);
        
        lower_bound(i) = find(cdf> 0.01,1, 'first');
        higher_bound(i) = find(cdf> 0.99,1, 'first');
        max_proj_rescaled(:,:,i)=mat2gray(max_proj{i},[lower_bound(i) higher_bound(i)]);
        
    end
    %          figure(handles.ax_main_image);
    imshow(max_proj_rescaled)
    
    BW_mask = roipoly;
    
    handles.stack.frame(handles.image_displayed).number_cells = handles.stack.frame(handles.image_displayed).number_cells + 1;
    handles.stack.frame(handles.image_displayed).cell(handles.stack.frame(handles.image_displayed).number_cells).mask  = BW_mask;
    
    if handles.stack.live_dot_finding
        [handles, dots] = dotCounting( handles, imdata, BW_mask);
        handles.stack.frame(handles.image_displayed).cell(handles.stack.frame(handles.image_displayed).number_cells).dots = dots;
    end
    close(h);

    guidata(hObject, handles)
    
    
    
end


% --------------------------------------------------------------------
function menu_removeCell_Callback(hObject, eventdata, handles)
% hObject    handle to menu_removeCell (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~isempty(handles.stack.frame(handles.image_displayed).cell)
    h= figure;
    num_channels = handles.num_channel;
    imdata = handles.current_image_stack;
    % Show max proj of current stack (first 3 channels)
    max_proj={};max_proj_rescaled=zeros(size(imdata{1},1),size(imdata{1},2),min(num_channels,3));
    for i=1:min(num_channels,3)
        max_proj{i} = max(imdata{i},[],3);
        
        % Automatic rescaling to display RGB
        bincounts = histc(max_proj{i}(:), [0.5:1:65535.5]);
        cdf=cumsum(bincounts)/size(max_proj{i},1)/size(max_proj{i},2);
        
        lower_bound(i) = find(cdf> 0.01,1, 'first');
        higher_bound(i) = find(cdf> 0.99,1, 'first');
        max_proj_rescaled(:,:,i)=mat2gray(max_proj{i},[lower_bound(i) higher_bound(i)]);
        
    end
    %          figure(handles.ax_main_image);
    imshow(max_proj_rescaled)
    
    [xClick, yClick, W]=ginput(1);
    xClick = round(xClick);
    yClick = round(yClick);
    click_image = zeros(size(max_proj_rescaled)) ;
    click_image(yClick,xClick) = 1;
    cellToRemove = 0;
    for i=1:handles.stack.frame(handles.image_displayed).number_cells
        match = find(logical(click_image).*handles.stack.frame(handles.image_displayed).cell(i).mask,1,'first');
        if match
            cellToRemove = i;
        end
    end
    if cellToRemove
        indexToKeep = setxor(1:handles.stack.frame(handles.image_displayed).number_cells,cellToRemove);
        handles.stack.frame(handles.image_displayed).number_cells = handles.stack.frame(handles.image_displayed).number_cells -1;
        handles.stack.frame(handles.image_displayed).cell = handles.stack.frame(handles.image_displayed).cell(indexToKeep);
    end
    close(h);
end
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
        handles.stack.image_path_number = [];
        handles.stack.frame.cell = struct;
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
% ---------------------CZI images-----------------------------------------------
function menu_addImages_Callback(hObject, eventdata, handles) %
% hObject    handle to menu_addImages (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[FileName,PathName] = uigetfile('*.czi','MultiSelect','on');

if ~isa(FileName, 'double') %Dialog not closed
    
    if isa(FileName, 'cell') % More than one file
       
        for i=1:length(FileName)
            handles.stack.image_path_cell{length( handles.stack.image_path_cell) + 1} = [PathName FileName{i}];
            handles.stack.frame(handles.stack.number_images + i).number_cells = 0;
        end
                
        handles.stack.number_images = handles.stack.number_images + length(FileName);
        
        
    else
        handles.stack.number_images = handles.stack.number_images + 1;
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
    prompt = ['Detected ' num2str(number_position) ' positions with ' num2str(number_channel) ' channels' newline 'Add positions (e.g. 1-5):'];
    dlg_title = 'Select positions';
    answer = inputdlg(prompt,dlg_title);
    if ~isempty(answer)
        str_cell = split(answer,'-');
        pos_start = str2num(str_cell{1});
        pos_end = str2num(str_cell{2});
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

function [imdata, num_channels] = tif_open(basename_path,actual_position)
files = dir([basename_path '*s' num2str(actual_position) '.tif']);
num_channels=length(files);
imdata={};
for i=1:num_channels
    imdata{i} = loadStack([files(i).folder filesep files(i).name]);
end


% --------------------------------------------------------------------
function handles = update_display(hObject,handles)
% TO DO: show cell outline
% TO DO: select cell to zoom on, keep focus on
if handles.image_displayed == 0 && handles.stack.number_images >= 1
    handles.image_displayed = 1;
end

if handles.image_displayed
    set(handles.txt_image_counter,'String', [num2str(handles.image_displayed) '/' num2str(handles.stack.number_images) ] );
    set(handles.slider_select_image, 'Value', handles.image_displayed);
    set(handles.slider_select_image, 'Min', 1);
    set(handles.slider_select_image, 'Max', handles.stack.number_images);
    if handles.stack.number_images >1
        set(handles.slider_select_image, 'SliderStep', [1/(handles.stack.number_images-1) , 1/(handles.stack.number_images-1) ]);
    else
        set(handles.slider_select_image, 'SliderStep', [0 , 0 ]);
    end
    if isempty(handles.current_image_stack)
        % Load the file stack 
        
        [handles.current_image_stack, handles.num_channel] = load_position(handles, handles.image_displayed);
            

    end
    num_channels = handles.num_channel;
    imdata = handles.current_image_stack;
    % Show max proj of current stack (first 3 channels)
    max_proj={};max_proj_rescaled=zeros(size(imdata{1},1),size(imdata{1},2),min(num_channels,3));
    for i=1:min(num_channels,3)
        max_proj{i} = max(imdata{i},[],3);
        
        % Automatic rescaling to display RGB
        bincounts = histc(max_proj{i}(:), [0.5:1:65535.5]);
        cdf=cumsum(bincounts)/size(max_proj{i},1)/size(max_proj{i},2);
        
        lower_bound(i) = find(cdf> 0.01,1, 'first');
        higher_bound(i) = find(cdf> 0.99,1, 'first');
        max_proj_rescaled(:,:,i)=mat2gray(max_proj{i},[lower_bound(i) higher_bound(i)]);
        
    end
    %          figure(handles.ax_main_image);
    imshow(max_proj_rescaled)
    
    if strcmp(get(handles.menu_stackViewer,'checked'),'on') %also update stackViewer
        handles = update_stackViewer_display(hObject,handles);
    end
     
end


% --------------------------------------------------------------------
function [handles, dots] = dotCounting( handles, imdata, cell_mask)
set(gcf,'pointer','watch');

BW_mask = cell_mask;
dots = struct;
for k = 1:handles.num_channel
    %             if k ~= nuc_channel_index
    im_stack = imdata{k};
    seg_im=uint16(BW_mask);

    masked_im_stack = bsxfun(@times, im_stack, seg_im);

    % Background substract (substract cell background)
    bg_level = median(nonzeros(masked_im_stack(:)));
    im_stack = max(0,im_stack - bg_level);
    
    
%     % 3D gaussian filtering
    SIGMA_GAUSS_FILT = 0.5;
    gauss_filt_imstack = imgaussfilt3(im_stack,SIGMA_GAUSS_FILT); %0*im_stack;
    pot_dots_im = imregionalmax(gauss_filt_imstack);
    log_filt_imstack = im_stack;
    for i =1:size(im_stack,3)
        log_filt_imstack(:,:,i) = logMask(im_stack(:,:,i));
        
    end
    

    
    
    masked_im_stack = bsxfun(@times, log_filt_imstack, seg_im);
    %some attempt at auto thresholding
    if handles.stack.autothreshold == 1
        frac_thres = graythresh(mat2gray(nonzeros(masked_im_stack)));
        thresh = frac_thres*(max(nonzeros(masked_im_stack))-min(nonzeros(masked_im_stack)))+min(nonzeros(masked_im_stack));
        handles.threshold.(['channel' num2str(k)]) = thresh;
        if strcmp(get(handles.menu_threshold,'Checked'), 'on') % Display threshold
            set(handles.threshold.(['edit_channel' num2str(k)]),'String',num2str(handles.threshold.(['channel' num2str(k)])));
        end
    else
        thresh = handles.threshold.(['channel' num2str(k)]);
        %TO DO: validate that we have a set threshold before using it
    end
    
    %                 thresh = 6000;
    thresh_stack = masked_im_stack  ;
    thresh_stack(thresh_stack < thresh) = 0;
    % alternative: take local maximum but don't threshold on area
    thresh_stack = double(thresh_stack) .* pot_dots_im;    
    
    conn_list = bwconncomp(thresh_stack,26);
    stats = regionprops(conn_list);
    idx_to_keep = [stats.Area] >= 1 ;
    good_dots_stats = stats(idx_to_keep);
    centroid_to_plot = zeros(length(good_dots_stats),3);
%     
%     for p = 1:length(good_dots_stats)
%         centroid_to_plot(p,:) = good_dots_stats(p).Centroid;
%     end
    
    dots(k).properties = good_dots_stats;
    dots(k).counts = length(good_dots_stats);
    
end
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


function lapFrame = logMask(im)   %**check that have the right logMask function

k = [-4 -1  0 -1 -4;...
     -1  2  3  2 -1;...
      0  3  4  3  0;...
     -1  2  3  2 -1;...
     -4 -1  0 -1 -4];

lapFrame = imfilter(im,k,'repl');


% --------------------------------------------------------------------
function menu_threshold_Callback(hObject, eventdata, handles)
% hObject    handle to menu_threshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


if strcmp(handles.menu_threshold.Checked,'on')
    handles.menu_threshold.Checked = 'off';
    menu_threshold_CloseRequest_Fcn(handles.threshold.window, 1, hObject);
else
    handles.menu_threshold.Checked = 'on';
    window_position = getpixelposition(handles.figure1);
    handles.threshold.window = figure('Name', 'Threshold settings', 'Menubar', 'none', 'Position', [window_position(1), window_position(2), 261, 271], 'CloseReq', {@menu_threshold_CloseRequest_Fcn, hObject}); % 'Keypress', {@onkeypress, hObject}'ResizeFcn', {@movie_slider_ResizeFcn, hObject});
    
    handles.threshold.chk_auto_threshold = uicontrol(gcf, 'Style', 'checkbox', 'Position', [31, 218, 95, 19], 'Value', 1, 'String', 'Auto Threshold', 'Callback', {@chk_auto_threshold_Callback, hObject});
    handles.stack.autothreshold = 1;

    handles.threshold.chk_live_dot_finding = uicontrol(gcf, 'Style', 'checkbox', 'Position', [31, 250, 95, 19], 'Value', 1, 'String', 'Live dot finding', 'Callback', {@chk_live_dot_finding_Callback, hObject});
    handles.stack.live_dot_finding = 1;
    
    handles.threshold.button_rethreshold = uicontrol(gcf, 'Style', 'pushbutton', 'Position', [141, 210, 84, 35], 'String', 'Rethreshold all', 'Callback', {@button_rethreshold_Callback, hObject});
    
    for i=1:handles.num_channel
        handles.threshold.(['edit_channel' num2str(i)]) = uicontrol(gcf, 'Style', 'edit', 'Position', [26, 170-40*(i-1), 46, 34], 'String', '', 'Callback', {@change_threshold_Callback, hObject}, 'Tag', ['edit_channel' num2str(i)]);
        handles.threshold.(['channel' num2str(i)]) = '';
    end
    %     handles.list_of_open_objects(end + 1) = handles.stackViewer.window;
end
guidata(hObject, handles)

function chk_auto_threshold_Callback(parent,~ ,hObject)

handles = guidata(hObject);

if get(parent,'value')
    handles.stack.autothreshold = 1;
else
    handles.stack.autothreshold = 0;
end

guidata(hObject, handles)

function chk_live_dot_finding_Callback(parent,~ ,hObject)

handles = guidata(hObject);

if get(parent,'value')
    handles.stack.live_dot_finding = 1;
else
    handles.stack.live_dot_finding = 0;
end

guidata(hObject, handles)

function button_rethreshold_Callback(parent,~ ,hObject)

handles = guidata(hObject);
h = waitbar(0);
for fov=1:handles.stack.number_images
   if  ~isempty(handles.stack.frame(fov).cell)
       for cellNo=1:handles.stack.frame(fov).number_cells
           handles = rethreshold_cell(handles, fov, cellNo);
           % TO DO: threshold all cells in fov at same time
       end
   end
   waitbar(fov/handles.stack.number_images,h);
end
close(h);
guidata(hObject, handles)

function change_threshold_Callback(parent,~ ,hObject)
handles = guidata(hObject);
parent_tag = get(parent,'Tag');
channel_to_change = parent_tag(end);
new_threshold = str2num(get(parent,'String'));

if ~handles.stack.autothreshold && ~isempty(new_threshold) %input validation
    handles.threshold.(['channel' channel_to_change]) = new_threshold;
    
else %put it back the way it was
   set(parent,'String',num2str(handles.threshold.(['channel' channel_to_change]))); 
end
% Rethreshold last cell
if handles.stack.frame(handles.image_displayed).number_cells && handles.stack.live_dot_finding
    handles = rethreshold_cell(handles, handles.image_displayed, length(handles.stack.frame(handles.image_displayed).cell));
    if strcmp(get(handles.menu_stackViewer,'checked'),'on') %also update display
        handles = update_stackViewer_display(hObject,handles);
    end
end
guidata(hObject, handles)


function menu_threshold_CloseRequest_Fcn(parent,~ ,hObject)
delete(parent);
handles = guidata(hObject);
handles.menu_threshold.Checked = 'off';
guidata(hObject, handles)


function handles = rethreshold_cell(handles, frameNo, cellNo)

 BW_mask = handles.stack.frame(frameNo).cell(cellNo).mask ; 
 if frameNo ~= handles.image_displayed %need to load the file
     imdata = load_position(handles,frameNo);
 else
     imdata = handles.current_image_stack;
 end
 [handles,dots] = dotCounting( handles, imdata, BW_mask);

  handles.stack.frame(frameNo).cell(cellNo).dots = dots;


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
     if strcmp(get(handles.menu_threshold,'Checked'),'on')
         menu_threshold_CloseRequest_Fcn(handles.threshold.window, 1, hObject);
     end
     

     delete(hObject);
   case 'No'
     return
 end







