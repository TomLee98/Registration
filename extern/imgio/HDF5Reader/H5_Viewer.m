function varargout = H5_Viewer2(varargin)
% H5_VIEWER2 MATLAB code for H5_Viewer2.fig
%      H5_VIEWER2, by itself, creates a new H5_VIEWER2 or raises the existing
%      singleton*.
%
%      H = H5_VIEWER2 returns the handle to a new H5_VIEWER2 or the handle to
%      the existing singleton*.
%
%      H5_VIEWER2('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in H5_VIEWER2.M with the given input arguments.
%
%      H5_VIEWER2('Property','Value',...) creates a new H5_VIEWER2 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before H5_Viewer2_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to H5_Viewer2_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help H5_Viewer2

% Last Modified by GUIDE v2.5 03-May-2017 11:29:03

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @H5_Viewer2_OpeningFcn, ...
                   'gui_OutputFcn',  @H5_Viewer2_OutputFcn, ...
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


% --- Executes just before H5_Viewer2 is made visible.
function H5_Viewer2_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to H5_Viewer2 (see VARARGIN)

% Choose default command line output for H5_Viewer2
handles.output = hObject;
[filename,pathname]  = uigetfile({'*.h5'});  
fname=fullfile(pathname,filename);
handles.Title.String=filename;
[piezo,volumes,lasers,t,img_idx,res]=LoadImgProperties(fname);
handles.Img=initialize_imgs(volumes,lasers,res);
[handles.Img,handles.T_vols]=get_all_volumes(fname,lasers,piezo,img_idx,res,volumes,t,handles.Img);
%handles.Img=flip(permute(handles.Img,[2,1,3,4,5]),1);
handles.Img_full=handles.Img;
handles.which_lasers=fliplr(find(fliplr(any(lasers,1))));
%handles.Img=double(handles.Img);
%handles.rgb=zeros([res(1),res(2),size(handles.Img,3),size(handles.Img,4),3]);
handles.TimeSlider.Min=1;
handles.TimeSlider.Max=size(handles.Img,4);
handles.TimeSlider.SliderStep=[1/size(handles.Img,4),3/size(handles.Img,4)];



handles.Z=round(size(handles.Img,3)/2);
handles.t=1;
handles.TimeSlider.Value=1;
handles.TimeText.String=sprintf('Frame %d/%d, t = %0.0f sec',handles.t,...
    size(handles.Img,4),handles.T_vols(handles.t));

handles.tplot=plot(handles.TimeAxis,handles.T_vols(handles.t)*ones(10,1),linspace(0,1,10));
handles.TimeAxis.XLim=[0,max(handles.T_vols)];

%img=double(handles.Img);
%handles.rgb(:,:,:,:,handles.which_lasers)=img;
%% Compute max intenstity projections
handles=max_inten_proj(handles);
%% scale image values
handles = scale_img_values(handles);
color={'Red','Green'};
for ii=1:length(color)
    
    if any(handles.which_lasers==ii)
        h=findobj('Tag',[color{ii},'Slider']);
            h.Min=0;
            which_idx=find(handles.which_lasers==ii);
            maxvalue=double(max(reshape(handles.Img(:,:,:,:,which_idx),[],1)));
            h.Max=log10(100*maxvalue);
    else
        h=findobj('Tag',[color{ii},'Slider']);
        h.Enable='off';
        h.Visible='off';
        hbutton=findobj('Tag',['Autoscale',color{ii}]);
        hbutton.Enable='off';
        hbutton.Visible='off';
    end
end

%handles.scale_rgb=handles.scale_rgb;
%handles.scale_rgb_slice=handles.scale_rgb_slice;
%handles.scale_rgb_proj=handles.scale_rgb_proj;
%image(handles.ax_slice,handles.Img(:,:,Z,t);
%image(handles.axes_max_proj,handles.Img_max_proj_disp(:,:,:,t))

frame_slice = zeros(size(handles.scale_rgb_slice));
frame_slice(:,:,handles.which_lasers)=double(squeeze(handles.Img(:,:,handles.Z,handles.t,:)));

frame_proj = zeros(size(handles.scale_rgb_proj));
frame_proj(:,:,handles.which_lasers) = double(squeeze(handles.Img_max_proj(:,:,handles.t,:)));
handles.crop=0;
handles.slice_img=image(handles.ax_slice,frame_slice./double(handles.scale_rgb_slice));
handles.ax_slice.XTick=[];
handles.ax_slice.YTick=[];
handles.proj_img=image(handles.axes_max_proj,frame_proj./double(handles.scale_rgb_proj));
handles.axes_max_proj.XTick=[];
handles.axes_max_proj.YTick=[];
handles=draw_z_slice(handles);

axis(handles.ax_slice,'equal');
axis(handles.axes_max_proj,'equal');




%handles=update_image(handles,Z,t);
% Update handles structure
set(gcf,'windowscrollWheelFcn',{@ZSlide,handles});guidata(hObject, handles);

% UIWAIT makes H5_Viewer2 wait for user response (see UIRESUME)
% uiwait(handles.figure1);

function handles=draw_z_slice(handles)
hold(handles.axes_max_proj,'on');
yrng=handles.axes_max_proj.YLim;
xrng=handles.axes_max_proj.XLim;
y1 = ones(100,1)*(size(handles.Img_max_proj,1)-handles.Z);
%y1=(y1-yrng(1))*size(handles.Img_max_proj,1)/(diff(yrng));
x1 = linspace(1,size(handles.Img,2),100);
%x1 = (x1-xrng(1))*size(handles.Img_max_proj,2)/diff(xrng);
if isfield(handles,'yz_line')
    delete(handles.yz_line);
end
handles.yz_line=plot(handles.axes_max_proj,x1,y1,'w');

y2 = linspace(1,size(handles.Img,1),100);
%y2=(y2-yrng(1))*size(handles.Img_max_proj,1)/(diff(yrng));
x2 =  ones(100,1)*(size(handles.Img_max_proj,2) - handles.Z);
%x2 = (x2-xrng(1))*size(handles.Img_max_proj,2)/diff(xrng);
if isfield(handles,'xz_line')
    delete(handles.xz_line);
end
handles.xz_line=plot(handles.axes_max_proj,x2,y2,'w');
hold(handles.axes_max_proj,'off');
zsl=size(handles.Img,3);
if zsl>1
    handles.ZSliceTxt.String=sprintf('Z Slice# %d / %d',handles.Z, zsl);
else
    handles.ZSliceTxt.String='2D Image';
end



function ZSlide(hObject, eventdata, handles)
UPDN = eventdata.VerticalScrollCount;
Z=handles.Z;
zsl=size(handles.Img,3);
Z = Z - UPDN;
        if (Z < 1)
            Z = 1;
        elseif (Z > zsl)
            Z = zsl;
        end

handles.Z=Z;

handles=update_image(handles,handles.Z,handles.t);


set(gcf,'windowscrollWheelFcn',{@ZSlide,handles});
guidata(hObject, handles);

1;

% --- Outputs from this function are returned to the command line.
function handles=update_image(handles,Z,t)
frame_slice = zeros(size(handles.scale_rgb_slice));
frame_slice(:,:,handles.which_lasers)=double(squeeze(handles.Img(:,:,Z,t,:)));

frame_proj = zeros(size(handles.scale_rgb_proj));
frame_proj(:,:,handles.which_lasers) = double(squeeze(handles.Img_max_proj(:,:,t,:)));

handles.slice_img.CData=frame_slice./handles.scale_rgb_slice;

handles.proj_img.CData=frame_proj./handles.scale_rgb_proj;
handles.ax_slice.XLim=[1,size(handles.Img,2)];
handles.ax_slice.YLim=[1,size(handles.Img,1)];


handles.axes_max_proj.XLim=[1,size(handles.Img_max_proj,2)];
handles.axes_max_proj.YLim=[1,size(handles.Img_max_proj,1)];
axis(handles.ax_slice,'equal');
axis(handles.axes_max_proj,'equal');
handles=draw_z_slice(handles);







function varargout = H5_Viewer2_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on slider movement.
function TimeSlider_Callback(hObject, eventdata, handles)
% hObject    handle to TimeSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
handles.t=round(hObject.Value);
handles.TimeText.String=sprintf('Frame %d/%d, t = %0.0f sec',handles.t,...
    size(handles.Img,4),handles.T_vols(handles.t));
delete(handles.tplot)
handles.tplot=plot(handles.TimeAxis,handles.T_vols(handles.t)*ones(10,1),linspace(0,1,10));
handles.TimeAxis.XLim=[0,max(handles.T_vols)];
handles.TimeAxis.YTick=[];
handles.TimeAxis.XTick=sort([handles.T_vols(handles.t),handles.TimeAxis.XTick]);
handles=update_image(handles,handles.Z,handles.t);
    set(gcf,'windowscrollWheelFcn',{@ZSlide,handles});guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function TimeSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TimeSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function RedSlider_Callback(hObject, eventdata, handles)
% hObject    handle to RedSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
    handles.scale_rgb(1)=10^hObject.Value;
    handles.scale_rgb_slice(:,:,1) =handles.scale_rgb(1);
    handles.scale_rgb_proj(:,:,1) =handles.scale_rgb(1);
    handles=update_image(handles,handles.Z,handles.t);
    set(gcf,'windowscrollWheelFcn',{@ZSlide,handles});guidata(hObject, handles);
    

% --- Executes during object creation, after setting all properties.
function RedSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to RedSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function GreenSlider_Callback(hObject, eventdata, handles)
% hObject    handle to GreenSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
    handles.scale_rgb(2)=10^hObject.Value;
    handles.scale_rgb_slice(:,:,2) = handles.scale_rgb(2);
    handles.scale_rgb_proj(:,:,2) = handles.scale_rgb(2);
    handles=update_image(handles,handles.Z,handles.t);
    set(gcf,'windowscrollWheelFcn',{@ZSlide,handles});
    guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function GreenSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to GreenSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in AutoscaleRed.
function AutoscaleRed_Callback(hObject, eventdata, handles)
% hObject    handle to AutoscaleRed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    idx=1;
    which_idx=find(handles.which_lasers==idx);
    handles.scale_rgb(idx)=max(reshape(handles.Img(:,:,:,handles.t,which_idx),[],1));
    handles.scale_rgb_slice(:,:,idx) = handles.scale_rgb(idx);
    handles.scale_rgb_proj(:,:,idx) = handles.scale_rgb(idx);
    handles=update_image(handles,handles.Z,handles.t);
    %handles.rgb_disp=handles.rgb(:,:,:,:,handles.which_lasers(ii))/handles.scale_rgb(ii);
    %handles.Img_max_proj_disp(:,:,ii,:)=handles.Img_max_proj(:,:,ii,:)/handles.scale_rgb(ii);
    set(gcf,'windowscrollWheelFcn',{@ZSlide,handles});
    guidata(hObject, handles);

% --- Executes on button press in AutoscaleGreen.
function AutoscaleGreen_Callback(hObject, eventdata, handles)
% hObject    handle to AutoscaleGreen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    idx=2;
    which_idx=find(handles.which_lasers==idx);
    handles.scale_rgb(idx)=max(reshape(handles.Img(:,:,:,handles.t,which_idx),[],1));
    handles.scale_rgb(idx)=max(reshape(handles.Img(:,:,:,handles.t,which_idx),[],1));
    handles.scale_rgb_slice(:,:,idx) = handles.scale_rgb(idx);
    handles.scale_rgb_proj(:,:,idx) = handles.scale_rgb(idx);
    handles.GreenSlider.Value = log10(handles.scale_rgb(idx));
    handles=update_image(handles,handles.Z,handles.t);
    %handles.rgb_disp=handles.rgb(:,:,:,:,handles.which_lasers(ii))/handles.scale_rgb(ii);
    %handles.Img_max_proj_disp(:,:,ii,:)=handles.Img_max_proj(:,:,ii,:)/handles.scale_rgb(ii);
    set(gcf,'windowscrollWheelFcn',{@ZSlide,handles});
    guidata(hObject, handles);

function handles=max_inten_proj(handles)
    img_size=size(handles.Img);
    handles.img_yz = flip(squeeze(max(handles.Img,[],2)),2);
    handles.img_xz = flip(permute(squeeze(max(handles.Img,[],1)),[2,1,3,4]),1);
    handles.img_xy = squeeze(max(handles.Img,[],3));
    handles.Img_max_proj = ones(img_size(1)+img_size(3)+1,...
        img_size(2)+img_size(3)+1,...
        img_size(4),size(handles.Img,5),'uint16') ; %3rd dimension is time
    handles.TimeAxis.YTick=[];
    handles.Img_max_proj(1:img_size(1),1:img_size(2),:,:) = handles.img_xy;
    handles.Img_max_proj(img_size(1)+2:end, 1:img_size(2),:,:) = handles.img_xz;
    handles.Img_max_proj(1:img_size(1), img_size(2)+2:end,:,:) = handles.img_yz;
    
    
    function handles = reshape_max_inten(handles)
        img_size=size(handles.Img);
        handles.Img_max_proj = ones(img_size(1)+img_size(3)+1,...
        img_size(2)+img_size(3)+1,...
        img_size(4),size(handles.Img,5),'uint16') ; %3rd dimension is time
        if handles.crop
            rect=handles.rect;
            handles.Img_max_proj(1:img_size(1),1:img_size(2),:,:) = ...
                handles.img_xy(rect(1):rect(1)+rect(3),rect(2):rect(2)+rect(4),:,:);
            handles.Img_max_proj(img_size(1)+2:end, 1:img_size(2),:,:) = ...
                handles.img_xz(:,rect(2):rect(2)+rect(4),:,:);
            handles.Img_max_proj(1:img_size(1), img_size(2)+2:end,:,:) = ...
                handles.img_yz(rect(1):rect(1)+rect(3),:,:,:);            
        else
            handles.Img_max_proj(1:img_size(1),1:img_size(2),:,:) = handles.img_xy;
            handles.Img_max_proj(img_size(1)+2:end, 1:img_size(2),:,:) = handles.img_xz;
            handles.Img_max_proj(1:img_size(1), img_size(2)+2:end,:,:) = handles.img_yz;
        end
    %handles.Img_max_proj_disp = handles.Img_max_proj;
    %proj_size=size(handles.Img_max_proj);
function handles = scale_img_values(handles)
    img_size=size(handles.Img);
    proj_size=size(handles.Img_max_proj);
    handles.scale_rgb=ones(3,1);
    handles.scale_rgb_slice=ones(img_size(1),img_size(2),3);
    handles.scale_rgb_proj = ones(proj_size(1),proj_size(2),3);
    %handles.rgb_disp=handles.rgb;
    for ii=1:length(handles.which_lasers)    
        handles.scale_rgb(handles.which_lasers(ii))=max(reshape(handles.Img(:,:,:,handles.Z,ii),[],1));
        handles.scale_rgb_slice(:,:,handles.which_lasers(ii)) = handles.scale_rgb(handles.which_lasers(ii));
        handles.scale_rgb_proj(:,:,handles.which_lasers(ii)) = handles.scale_rgb(handles.which_lasers(ii));
        %handles.rgb_disp=handles.rgb(:,:,:,:,handles.which_lasers(ii))/handles.scale_rgb(ii);
        %handles.Img_max_proj_disp(:,:,ii,:)=handles.Img_max_proj(:,:,ii,:)/handles.scale_rgb(ii);
    end
function handles = reshape_scale_img_values(handles)
    img_size=size(handles.Img);
    proj_size=size(handles.Img_max_proj);
    %handles.scale_rgb=ones(3,1);
    handles.scale_rgb_slice=ones(img_size(1),img_size(2),3);
    handles.scale_rgb_proj = ones(proj_size(1),proj_size(2),3);
    %handles.rgb_disp=handles.rgb;
    for ii=1:length(handles.which_lasers)    
        %handles.scale_rgb(handles.which_lasers(ii))=max(reshape(handles.Img(:,:,:,handles.Z,ii),[],1));
        handles.scale_rgb_slice(:,:,handles.which_lasers(ii)) = handles.scale_rgb(handles.which_lasers(ii));
        handles.scale_rgb_proj(:,:,handles.which_lasers(ii)) = handles.scale_rgb(handles.which_lasers(ii));
        %handles.rgb_disp=handles.rgb(:,:,:,:,handles.which_lasers(ii))/handles.scale_rgb(ii);
        %handles.Img_max_proj_disp(:,:,ii,:)=handles.Img_max_proj(:,:,ii,:)/handles.scale_rgb(ii);
    end
% --- Executes during object creation, after setting all properties.
function ax_slice_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ax_slice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate ax_slice


% --------------------------------------------------------------------
function Untitled_1_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in Crop.
function Crop_Callback(hObject, eventdata, handles)
% hObject    handle to Crop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~handles.crop
    handles.slice_img.Visible='off';
   
    rect=getrect(handles.axes_max_proj);
    handles.slice_img.Visible='on';
    if all(rect(3:4)>0) 
        if rect(1)<1
           rect(1)=1;
        end
        if rect(2)<1
            rect(2)=1;
        end
        if rect(1)+rect(3)>size(handles.Img,2)
            rect(3)=size(handles.Img,1)-rect(1);
        end
        if rect(2)+rect(4)>size(handles.Img,1)
            rect(4)=size(handles.Img,2)-rect(2);
        end
        rect=rect([2,1,4,3]);
        rect=round(rect);
        handles.rect=rect;
        handles.Img=handles.Img(rect(1):rect(1)+rect(3),rect(2):rect(2)+rect(4),...
            :,:,:);
        handles.crop=1;
        handles = max_inten_proj(handles);
        handles = reshape_scale_img_values(handles);
        handles = update_image(handles,handles.Z,handles.t);
        
        handles.Crop.String='Restore';
        set(gcf,'windowscrollWheelFcn',{@ZSlide,handles});
        guidata(hObject, handles);
        
    end
else
    handles.Img=handles.Img_full;
    handles.crop=0;
    handles = max_inten_proj(handles);
    handles = reshape_scale_img_values(handles);
    handles = update_image(handles,handles.Z,handles.t);
    
    handles.Crop.String='Zoom';
    set(gcf,'windowscrollWheelFcn',{@ZSlide,handles});
    guidata(hObject, handles);    
end
1;
