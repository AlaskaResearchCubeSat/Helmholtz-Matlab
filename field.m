function varargout = field(varargin)
% FIELD Run Helmholtz cage GUI
%      FIELD, by itself, creates opens a new GUI if none exists or raises
%      the current window.
%
% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @field_OpeningFcn, ...
                   'gui_OutputFcn',  @field_OutputFcn, ...
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


% --- Executes just before field is made visible.
function field_OpeningFcn(hObject, ~, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to field (see VARARGIN)

% Choose default command line output for field
handles.output = hObject;

%check if GUI is already open
if ~isfield(handles,'Hobj')
    %initialize the helmholtz cage hardware
    handles.Hobj=cage_control();

    if(size(varargin,2)~=0)
        %get transform matrix
        T=varargin{1};
        %if passed an empty matrix then calibrate
        if(size(T)~=[4 4])
            %run calibration
            handles.Hobj.calibrate();
        else
            %set calibration
            handles.Hobj.T=T;
        end
    end

    %set field to zero
    handles.Hobj.Bs=[0 0 0];


    %set initial field plot
    set_graph(handles,[0 0 0]);
    %set initial view
    view(handles.field_graph,30,40);

    handles.simTimer = timer('TimerFcn',{@update,handles.Hobj},'Period',1,'ExecutionMode','fixedRate');

    handles.fieldFileData=[];

    %TODO: maybe read saved paths from file?
    
    %path for calibration files
    handles.calPath='C:\helmholtz-matlab\';
    %path for field sequence files
    handles.fieldPath='Z:\Helmholtz-Cage\Helmholtz-Programs\';
end
    
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes field wait for user response (see UIRESUME)
% uiwait(handles.field);

% --- Outputs from this function are returned to the command line.
function varargout = field_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.Hobj.T;


% --- Executes on button press in cal.
function cal_Callback(hObject, eventdata, handles)
% hObject    handle to cal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

cerr=handles.Hobj.calibrate();
%reset field to input values
set_calc_field(handles);



function xfield_Callback(hObject, eventdata, handles)
% hObject    handle to xfield (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of xfield as text
%        str2double(get(hObject,'String')) returns contents of xfield as a double


% --- Executes during object creation, after setting all properties.
function xfield_CreateFcn(hObject, eventdata, handles)
% hObject    handle to xfield (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function yfield_Callback(hObject, eventdata, handles)
% hObject    handle to yfield (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of yfield as text
%        str2double(get(hObject,'String')) returns contents of yfield as a double


% --- Executes during object creation, after setting all properties.
function yfield_CreateFcn(hObject, eventdata, handles)
% hObject    handle to yfield (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function zfield_Callback(hObject, eventdata, handles)
% hObject    handle to zfield (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of zfield as text
%        str2double(get(hObject,'String')) returns contents of zfield as a double

% --- Executes during object creation, after setting all properties.
function zfield_CreateFcn(hObject, eventdata, handles)
% hObject    handle to zfield (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in cal_tst.
function cal_tst_Callback(hObject, eventdata, handles)
% hObject    handle to cal_tst (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[max_err,mean_err]=handles.Hobj.cal_test();
set_calc_field(handles);
if ~isnan(max_err)
    msgbox(sprintf('The calibration test is complete.\nMean error = %.4g gauss\nMax error = %.4g',mean_err,max_err),'Calibration Test Complete');
end

%set the location edit values
function set_loc(handles,lat,long,R)
%set edit values
set(handles.lat,'String',sprintf('%.2f',lat));
set(handles.long,'String',sprintf('%.2f',long));
set(handles.alt,'String',sprintf('%.2f',R));

%set spherical coordinate edits from field
function set_spherical(handles,B)
    %calculate magnitude
    r=sqrt(sum(B.^2));
    %calculate phi
    phi=180/pi*acos(B(3)/r);
    %calculate theta
    theta=180/pi*atan2(B(2),B(1));
    %if r is zero then theta will be NaN, fix this
    if(r==0)
        %set phi to zero
        phi=0;
    end
    %set edit values
    set(handles.phi,'String',sprintf('%.2f',phi));
    set(handles.theta,'String',sprintf('%.2f',theta));
    set(handles.mag,'String',sprintf('%.2f',r));

%set rectangular coordinate edits from field
function set_rec(handles,B)
    %this is easy, just shove'm in
    set(handles.xfield,'String',sprintf('%.2f',B(1)));
    set(handles.yfield,'String',sprintf('%.2f',B(2)));
    set(handles.zfield,'String',sprintf('%.2f',B(3)));  
    
%set field graph from magnetic field vector
function set_graph(handles,B)
    %clear old plot
    cla(handles.field_graph)
    %set hold
    hold(handles.field_graph,'on');
    %plot X component
    plot3(handles.field_graph,[0 B(1)],[0 0],[0 0],'r');
    %plot Y component
    plot3(handles.field_graph,[0 0],[0 B(2)],[0 0],'g');
    %plot Z component
    plot3(handles.field_graph,[0 0],[0 0],[0 B(3)],'b');
    %plot field vector
    plot3(handles.field_graph,[0 B(1)],[0 B(2)],[0 B(3)],'k');
    %unset hold
    hold(handles.field_graph,'off');
    %display legend
    legend(handles.field_graph,'X','Y','Z','Magnitude');
    legend(handles.field_graph,'Location','NorthEast');
    %get limit for graph
    lim=1.1*max([max(B) 20]);
    lim=min([lim 2000]);
    %set axis to a square centered at zero
    axis(handles.field_graph,[-lim lim -lim lim -lim lim]);
    %allow the graph to be rotated
    rotate3d(handles.field_graph,'on');
    %set axis labels
    xlabel(handles.field_graph,'X - axis Magnetic Field [miligauss]');
    ylabel(handles.field_graph,'Y - axis Magnetic Field [miligauss]');
    zlabel(handles.field_graph,'Z - axis Magnetic Field [miligauss]');

%set field from value in edit boxes
function B=set_calc_field(handles)
%empty field vector
B=zeros(3,1);
%get the spherical button
b=findobj('Tag','spherical');
if get(b,'Value') 
    %use spherical coordinates
    phi=str2double(get(handles.phi,'String'))*pi/180;
    if isnan(phi)
        phi=0;
        set(handles.phi,'String','0');
    end
    theta=str2double(get(handles.theta,'String'))*pi/180;
    if isnan(theta)
        theta=0;
        set(handles.theta,'String','0');
    end
    r=str2double(get(handles.mag,'String'));
    if isnan(r)
        r=0;
        set(handles.mag,'String','0');
    end
    %calculate field vector in rectangular coordinates
    B(1)=r*cos(theta)*sin(phi);
    B(2)=r*sin(theta)*sin(phi);
    B(3)=r*cos(phi);
    %set rectangular coordinate edits
    set_rec(handles,B);
    %set field graph
    set_graph(handles,B);
else
    %use rectangular coordinates
    B(1)=str2double(get(handles.xfield,'String'));
    if isnan(B(1))
        B(1)=0;
        set(handles.xfield,'String','0');
    end
    B(2)=str2double(get(handles.yfield,'String'));
    if isnan(B(2))
        B(2)=0;
        set(handles.yfield,'String','0');
    end
    B(3)=str2double(get(handles.zfield,'String'));
    if isnan(B(3))
        B(3)=0;
        set(handles.zfield,'String','0');
    end
    %set spherical coordinate edits
    set_spherical(handles,B);
    %set field graph
    set_graph(handles,B);
end
%convert to gauss
B=B*1e-3;
%make sure everything went ok
if isnan(B)
    fprintf(1,'Error: with field values\n');
else
    %set field in Helmholtz cage
    handles.Hobj.Bs=B;
end

% --- Executes on button press in set_field.
function set_field_Callback(hObject, eventdata, handles)
% hObject    handle to set_field (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(hObject,'Enable','off');
set_calc_field(handles);
set(hObject,'Enable','on');

% --- Executes when user attempts to close field.
function field_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to field (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

Hobj=handles.Hobj;

stop(handles.simTimer);
delete(handles.simTimer);

% Hint: delete(hObject) closes the figure
delete(hObject);
delete(Hobj);

% --- Executes during object creation, after setting all properties.
function plot_CreateFcn(hObject, eventdata, handles)
% hObject    handle to plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% --- Executes when field is resized.
function field_ResizeFcn(hObject, eventdata, handles)
% hObject    handle to field (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%TODO: should write something here

function phi_Callback(hObject, eventdata, handles)
% hObject    handle to phi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of phi as text
%        str2double(get(hObject,'String')) returns contents of phi as a double


% --- Executes during object creation, after setting all properties.
function phi_CreateFcn(hObject, eventdata, handles)
% hObject    handle to phi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function theta_Callback(hObject, eventdata, handles)
% hObject    handle to theta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of theta as text
%        str2double(get(hObject,'String')) returns contents of theta as a double

% --- Executes during object creation, after setting all properties.
function theta_CreateFcn(hObject, eventdata, handles)
% hObject    handle to theta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function mag_Callback(hObject, eventdata, handles)
% hObject    handle to mag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mag as text
%        str2double(get(hObject,'String')) returns contents of mag as a
%        double

% --- Executes during object creation, after setting all properties.
function mag_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%enable UI controls based on current state of checkboxes
function setCoordEn(handles)
sph={'phi','theta','mag'};
rec={'xfield','yfield','zfield'};
fctl={'rectangular','spherical','set_field'};
loc={'lat','long','alt','set_loc','cur_date'};
lctl={'loc2field'};
date={'date'};

orb_parm={'incTxt','inclination','incUnit','altTxt','altitude','altUnit',...
    'raanTxt','raan','raanUnit','simSpdTxt','sim_spd_s','sim_spd_e'};
orb_file={'file_name_disp','file'};
orb={'sim_run','orb_use_file','orb_use_parm','stepTxt','sim_step','stepUnit'};

%check the set location from field check box
if get(handles.orb2loc,'Value')
    %check if simulation is running
    if get(handles.sim_run,'Value')
        enctl={'sim_run'};
        disctl={rec{:},sph{:},fctl{:},date{:},loc{:},lctl{:},orb_file{:},orb_parm{:},orb{:}};
    else
        if get(handles.orb_use_parm,'Value')
            %use orbital paramiters
            enctl={orb{:},orb_parm{:}};
            disctl={rec{:},sph{:},fctl{:},date{:},loc{:},lctl{:},orb_file{:}};
        else
            %use file
            enctl={orb{:},orb_file{:}};
            disctl={rec{:},sph{:},fctl{:},date{:},loc{:},lctl{:},orb_parm{:}};
        end
    end
else
    %check the set field from location check box
    if get(handles.loc2field,'Value')
        %check the current date check box
        if get(handles.cur_date,'Value')
            enctl={loc{:},lctl{:}};
            disctl={rec{:},sph{:},fctl{:},date{:},orb{:},orb_file{:},orb_parm{:}};
        else
            enctl={loc{:},date{:},lctl{:}};
            disctl={rec{:},sph{:},fctl{:},orb{:}};    
        end
    else
        %check the spherical/rectangular radio button
        if get(handles.spherical,'Value') 
            enctl={sph{:},fctl{:},lctl{:}};
            disctl={rec{:},loc{:},date{:},orb{:},orb_file{:},orb_parm{:}};
        else
            enctl={rec{:},fctl{:},lctl{:}};
            disctl={sph{:},loc{:},date{:},orb{:},orb_file{:},orb_parm{:}};
        end
    end
end
%disable all controls in disctl
for k=disctl
    b=findobj('Tag',k{1});
    set(b,'Enable','off');
end
%enable all controls in enctl
for k=enctl
    b=findobj('Tag',k{1});
    set(b,'Enable','on');
end   

% --- Executes on button press in loc2field.
function loc2field_Callback(hObject, eventdata, handles)
% hObject    handle to loc2field (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of loc2field

setCoordEn(handles)

function lat_Callback(hObject, eventdata, handles)
% hObject    handle to lat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of lat as text
%        str2double(get(hObject,'String')) returns contents of lat as a double

% --- Executes during object creation, after setting all properties.
function lat_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function long_Callback(hObject, eventdata, handles)
% hObject    handle to long (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of long as text
%        str2double(get(hObject,'String')) returns contents of long as a double


% --- Executes during object creation, after setting all properties.
function long_CreateFcn(hObject, eventdata, handles)
% hObject    handle to long (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function alt_Callback(hObject, eventdata, handles)
% hObject    handle to alt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of alt as text
%        str2double(get(hObject,'String')) returns contents of alt as a double

% --- Executes during object creation, after setting all properties.
function alt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to alt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes when selected object is changed in ctls.
function ctls_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in ctls 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
tag=get(eventdata.NewValue,'Tag');
switch  tag% Get Tag of selected object.
    case 'spherical'
        setCoordEn(handles);
    case 'rectangular'
        setCoordEn(handles);
    otherwise
        %there is no match.
        printf('Tag "%s"\n',tag);
end

% --- Executes on button press in set_loc.
function set_loc_Callback(hObject, eventdata, handles)
% hObject    handle to set_loc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

lat=str2double(get(handles.lat,'String'));
long=str2double(get(handles.long,'String'));
alt=str2double(get(handles.alt,'String'));
date=str2double(get(handles.date,'String'));
%get magnetic field using magnetic field model
try
    B=magfd(date,2,alt,90-lat,long,eye(3))*1e4*1e3;%convert to miligauss
catch
    fprintf(2,'Error calculating magnetic field, using 2010 for date.\n')
    B=magfd(2010,2,alt,90-lat,long,eye(3))*1e4*1e3;%convert to miligauss
end
%set edit boxes to value from model
set_rec(handles,B);
set_spherical(handles,B);
set_graph(handles,B);
%convert to gauss
B=B*1e-3;
if isnan(B)
    fprintf(1,'Error: with field values\n');
else
    handles.Hobj.Bs=B;
end

function date_Callback(hObject, eventdata, handles)
% hObject    handle to date (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of date as text
%        str2double(get(hObject,'String')) returns contents of date as a double

% --- Executes during object creation, after setting all properties.
function date_CreateFcn(hObject, eventdata, handles)
% hObject    handle to date (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
set(hObject,'String',sprintf('%4.2f',now/365.24));

% --- Executes on button press in cur_date.
function cur_date_Callback(hObject, eventdata, handles)
% hObject    handle to cur_date (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cur_date
if(get(hObject,'Value'))
    set(handles.date,'String',sprintf('%4.2f',now/365.24));
end
setCoordEn(handles)

% --- Executes on button press in off_butt.
function off_butt_Callback(hObject, eventdata, handles)
% hObject    handle to off_butt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.Hobj.off();

% --- Executes on button press in orb2loc.
function orb2loc_Callback(hObject, eventdata, handles)
% hObject    handle to orb2loc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

setCoordEn(handles)
set(handles.sim_run,'Value',0);
%tell run button that value has changed
sim_run_Callback(handles.sim_run,eventdata,handles);

% --- Executes on slider movement.
function sim_spd_s_Callback(hObject, eventdata, handles)
% hObject    handle to sim_spd_s (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
val=get(hObject,'Value');
set(handles.sim_spd_e,'String',sprintf('%.2f',val));

% --- Executes during object creation, after setting all properties.
function sim_spd_s_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sim_spd_s (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

% --- Executes on button press in sim_run.
function sim_run_Callback(hObject, eventdata, handles)
% hObject    handle to sim_run (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if get(hObject,'Value')
    %read simulation step length
    step=str2double(get(handles.sim_step,'String'));
    %set timer period to simulation step length
    handles.simTimer.Period=step;
    %check if file or parameters are used
    if get(handles.orb_use_parm,'Value')
        %generate data from orbital parameters
        speed=str2double(get(handles.sim_spd_e,'String'));
        raan=str2double(get(handles.raan,'String'));
        alt=str2double(get(handles.altitude,'String'));
        inc=str2double(get(handles.inclination,'String'));
        %get date
        date=str2double(get(handles.date,'String'));
        %generate UserData structure
        handles.simTimer.UserData=struct('step',step*speed,'time',0,...
            'alt',alt,'raan',raan,'inc',inc','date',date,'handles',handles,'file',false);
        %check the set field from location check box
        set(handles.loc2field,'Value',true);
    else
        size(handles.fieldFileData)
        %check data to see if it is present
        if isempty(handles.fieldFileData)
            %no data clear checkbox and return
            set(hObject,'Value',0);
            %return
            return;
        end
        %generate UserData structure
        handles.simTimer.UserData=struct('idx',1,'handles',handles,'file',true,'B',handles.fieldFileData);
    end
    %disable calibration and calibration test
    set(handles.cal,'Enable','off');
    set(handles.cal_tst,'Enable','off');
    %start timer
    start(handles.simTimer);
else
    %stop timer
    stop(handles.simTimer);
    %enable calibration and calibration test
    set(handles.cal,'Enable','on');
    set(handles.cal_tst,'Enable','on');
end
%set enabled/disabled status of input controls
setCoordEn(handles);

function sim_spd_e_Callback(hObject, eventdata, handles)
% hObject    handle to sim_spd_e (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sim_spd_e as text
%        str2double(get(hObject,'String')) returns contents of sim_spd_e as
%        a double

%get value
val=get(hObject,'String');
%convert value to double
val=str2double(val);
%check if conversion worked
if ~isnan(val)
    %get minimum value
    min=get(handles.sim_spd_s,'Min');
    %set lower bound to min
    if val<min
        val=min;
        %set new edit value
        set(hObject,'String',sprintf('%.2f',val));
    end
    %get maximum value
    max=get(handles.sim_spd_s,'Max');
    %set upper bound to max
    if val>max
       val=max;
       %set new edit value
        set(hObject,'String',sprintf('%.2f',val));
    end 
    %set scroll bar value
    set(handles.sim_spd_s,'Value',val);
end

% --- Executes during object creation, after setting all properties.
function sim_spd_e_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sim_spd_e (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in on_butt.
function on_butt_Callback(hObject, eventdata, handles)
% hObject    handle to on_butt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.Hobj.on();

% --- Executes on button press in save_cal.
function save_cal_Callback(hObject, eventdata, handles)
% hObject    handle to save_cal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[name,path,type]=uiputfile('*.cal','Save Calibration','calibration.cal');
if not(isequal(name,0) && isequal(path,0) && isequal(type,0))
    handles.Hobj.saveCal([path name])
    %save calibration path
    handles.calPath=path;
end


% --- Executes on button press in load_cal.
function load_cal_Callback(hObject, eventdata, handles)
% hObject    handle to load_cal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[name,path]=uigetfile('*.cal','Load Calibration',[handles.calPath 'calibration.cal']);
if not(isequal(name,0) && isequal(path,0))
    
    try
        %load new calibration
        handles.Hobj.loadCal([path name]);
        %save calibration path
        handles.calPath=path;
    catch
        errordlg(['Error reading calibration file ''' path name '''. Pplease check to see if the file is in the proper format.'],'File Error','modal');
        return;
    end
    %set field based on new calibration
    set_calc_field(handles);
end


% --- Executes on button press in reset_view.
function reset_view_Callback(hObject, eventdata, handles)
% hObject    handle to reset_view (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Orient the view to azimuth = 30 and elevation = 40.
view(handles.field_graph,30,40);

% --- Executes on button press in show_cal.
function show_cal_Callback(hObject, eventdata, handles)
% hObject    handle to show_cal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
msgbox([sprintf('Calibration Matrix\n') sprintf('%7.2f %7.2f %8.2f %7.2f\n',handles.Hobj.T')],'Calibration Matrix');


% --- Executes on button press in file.
function file_Callback(hObject, eventdata, handles)
% hObject    handle to file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[name,path]=uigetfile('*','Load File',handles.fieldPath);
if not(isequal(name,0) && isequal(path,0))
    try
        %load data from file
        dat=load([path name]);
        %set path for field files
        handles.fieldPath=path;
    catch
        errordlg(['Error reading file ''' path name '''. please check to see if the file is in the proper format.'],'File Error','modal');
        return;
    end 
    %set data array
    handles.fieldFileData=dat;
    %set file name display
    set(handles.file_name_disp,'String',name);
    %update data structure
    guidata(hObject, handles);
end


% --- Executes when selected object is changed in orb_group.
function orb_group_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in orb_group 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)

%new=get(eventdata.NewValue,'String');
%old=get(eventdata.OldValue,'String');
%fprintf(1,'Old Value %s\nNew Value %s\n',old,new);
%get(eventdata.NewValue,'String')
setCoordEn(handles)

%switch new
%    case 'Run'
        



function inclination_Callback(hObject, eventdata, handles)
% hObject    handle to inclination (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of inclination as text
%        str2double(get(hObject,'String')) returns contents of inclination as a double


% --- Executes during object creation, after setting all properties.
function inclination_CreateFcn(hObject, eventdata, handles)
% hObject    handle to inclination (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function altitude_Callback(hObject, eventdata, handles)
% hObject    handle to altitude (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of altitude as text
%        str2double(get(hObject,'String')) returns contents of altitude as a double


% --- Executes during object creation, after setting all properties.
function altitude_CreateFcn(hObject, eventdata, handles)
% hObject    handle to altitude (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function raan_Callback(hObject, eventdata, handles)
% hObject    handle to raan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of raan as text
%        str2double(get(hObject,'String')) returns contents of raan as a double


% --- Executes during object creation, after setting all properties.
function raan_CreateFcn(hObject, eventdata, handles)
% hObject    handle to raan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function epoc_Callback(hObject, eventdata, handles)
% hObject    handle to epoc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of epoc as text
%        str2double(get(hObject,'String')) returns contents of epoc as a double


% --- Executes during object creation, after setting all properties.
function epoc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to epoc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
set(hObject,'String',sprintf('%4.2f',now/365.24));



function sim_step_Callback(hObject, eventdata, handles)
% hObject    handle to sim_step (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sim_step as text
%        str2double(get(hObject,'String')) returns contents of sim_step as a double


% --- Executes during object creation, after setting all properties.
function sim_step_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sim_step (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%update field callback called from orbit timer
function update(obj,event,Hobj)
    %get data structure from event user data
    dat=obj.UserData;
    if dat.file
        %advance to next index
        dat.idx=dat.idx+1;
        %check bounds
        if dat.idx>length(dat.B)
            %reset index
            dat.idx=1;
        end
        %write back structure into UserData
        obj.UserData=dat;
        %read flux vecotr from array
        FluxVector=dat.B(dat.idx,:);
    else
        %update simulation time
        dat.time=dat.step+dat.time;
        %write back structure into UserData
        obj.UserData=dat;
        %calculate latitude and longitude
        [lat,long]=calcLatLong(dat.time,dat.inc,dat.raan,dat.alt);
        % Radius of the Earth (m)[IUGG value of equatorial radius]; Orbit Altitude(m)
        Re=6.378137*10^6;
        % Position vector magnitude (km)
        R=dat.alt+(Re/1e3); 
        %calculate field vector based on location
        [FluxVector] = magfd(2012.5,2,R,90-lat,long,eye(3))*1e4;
    end
    %set field
    Hobj.Bs=FluxVector;
    
    %convert to miligauss for display
    FluxVector=FluxVector*1e3;
    %set spherical coordinat edits
    set_spherical(dat.handles,FluxVector);
    %set rectangular coordinat edits
    set_rec(dat.handles,FluxVector);
    %set smagnetic field graph
    set_graph(dat.handles,FluxVector);
    if ~dat.file
        %set location (latitude, longitude, altitude)
        set_loc(dat.handles,lat,long,R);
    end
    
 %calculate latitude and longitude from orbit parms
 %most of this code was lifted from Cranks thesis
 %time in secconds
 %inc and raan in deg
 %alt in km
function [lat,long]=calcLatLong(time,inc,raan,alt)
    % Gravitational constant (m^3/sec^2)
    GM = 3.986004418*10^14;
    % Radius of the Earth (m)[IUGG value of equatorial radius]; Orbit Altitude(m)
    Re = 6.378137*10^6;
    % Position vector magnitude (km)
    R =  alt+(Re/1e3);            
    %convert to radians
    inc=inc*pi/180;
    lambda_not=raan*pi/180;
    % Rate that lambda changes (rad/sec)
    lambda_dot = -4.17089*10^7 *(R)^(-7/2)*cos(inc); 
    % Orbital Rate (rad/sec)
    CapOmega_not=sqrt(GM/(R*1e3)^3);   
    % Orbital period (sec)
    %Tp = (2*pi)/CapOmega_not;  
    % Earth rotation rate (2 *pi in 24 hours*3600sec/hour)
    gamma_dot = (2*pi)/(24*3600);  
    %don't really know what these is but set to zero
    eta_not=0;% Initial S/C Position wrt ascending node (degrees); 
    gamma_not = 0; % Initial Magnetic axis position wrt Inertial
    % Tracks current spacecraft position wrt ascending node
    %etaO = eta_not + CapOmega_not*time;         

    
    % Calculates current magnetic reference angle
    %used to relate ECI & ECEF frames
    gamma = gamma_not + gamma_dot*time; 
    % Calculates current spacecraft position wrt the ascending node
    eta = eta_not + CapOmega_not*time;
    % Calculates current position of the longitude of the ascending node 
    lambda = lambda_not + lambda_dot*time;

    % Calcs required to convert eta, gamma, incline, lambda into latitude
    % and longitude
    ui=R*cos(eta)-R*sin(eta)*cos(inc);
    vi=R*cos(eta)+R*sin(eta)*cos(inc);
    we=R*sin(eta)*sin(inc);
    ue=cos(gamma+lambda)*ui+sin(gamma+lambda)*vi;
    ve=-sin(gamma+lambda)*ui+cos(gamma+lambda)*vi;
    lat = asin(we/R)*180/pi;
    long = (ve/abs(ve))*acos(ue/sqrt(ue^2+ve^2))*180/pi-45;
    
