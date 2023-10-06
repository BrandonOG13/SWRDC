function varargout = gui(varargin)
% GUI MATLAB code for gui.fig
%      GUI, by itself, creates a new GUI or raises the existing
%      singleton*.
%
%      H = GUI returns the handle to a new GUI or the handle to
%      the existing singleton*.
%
%      GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI.M with the given input arguments.
%
%      GUI('Property','Value',...) creates a new GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help gui

% Last Modified by GUIDE v2.5 21-Oct-2013 22:06:34

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @gui_OpeningFcn, ...
                   'gui_OutputFcn',  @gui_OutputFcn, ...
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


% --- Executes just before gui is made visible.
function gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to gui (see VARARGIN)

% Choose default command line output for gui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes when figure1 is resized.
function figure1_ResizeFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
%Record command window
r = randi([1 1000],1);
from=clock;
diary(['commandwindow\log',...
    num2str(from(1)),'_',num2str(from(2)),'_',num2str(from(3)),...
    '_',num2str(from(4)),'_',num2str(from(5)),'_',num2str(from(6)),'.txt'])

%Collect input data
[data] = collect_data(handles);

%If endPop.mat will be loaded, check if it exists, If not, return to GUI
if get(handles.flag3,'Value')==1
    try
        load endPop.mat
    catch
        disp('endPop.mat file is missing or is not placed in the correct directory.')
        diary off
        clear all
        return
    end
end

%Check if input data is correct, else return to GUI
if isstruct(data)
    
    %Determine if the code is to be executed in parallel and if so, setup the
    %number of cores
    if get(handles.parallel,'Value')==1
        numcores=get(handles.numcores,'String');
        numcores = str2num(numcores);
        delete(gcp('nocreate'));
        try % Check if matlabpool works
            parpool('local', numcores);
        catch err % If not, display the matlab error message
            fprintf('\n');
            fprintf('Please read the MATLAB error message regarding your \n');
            fprintf('attempt to run SWRDC in parallel. \n');
            fprintf('\n');
            rethrow(err)
        end
    end
    
    [P,I2]=ga_ms(data); %Execute the optimization
    
    if get(handles.parallel,'Value')==1
        delete(gcp('nocreate'));
    end
else
    diary off
    clear all
    return
end

%Stop recording command window
diary off
clear all



function iter_Callback(hObject, eventdata, handles)
% hObject    handle to iter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of iter as text
%        str2double(get(hObject,'String')) returns contents of iter as a double
handles.data.iter=hObject;
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function iter_CreateFcn(hObject, eventdata, handles)
% hObject    handle to iter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function belements_Callback(hObject, eventdata, handles)
% hObject    handle to belements (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of belements as text
%        str2double(get(hObject,'String')) returns contents of belements as a double


% --- Executes during object creation, after setting all properties.
function belements_CreateFcn(hObject, eventdata, handles)
% hObject    handle to belements (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function u_start_Callback(hObject, eventdata, handles)
% hObject    handle to u_start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of u_start as text
%        str2double(get(hObject,'String')) returns contents of u_start as a double


% --- Executes during object creation, after setting all properties.
function u_start_CreateFcn(hObject, eventdata, handles)
% hObject    handle to u_start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function lambda_start_Callback(hObject, eventdata, handles)
% hObject    handle to lambda_start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of lambda_start as text
%        str2double(get(hObject,'String')) returns contents of lambda_start as a double


% --- Executes during object creation, after setting all properties.
function lambda_start_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lambda_start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Jgen_Callback(hObject, eventdata, handles)
% hObject    handle to Jgen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Jgen as text
%        str2double(get(hObject,'String')) returns contents of Jgen as a double


% --- Executes during object creation, after setting all properties.
function Jgen_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Jgen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function tor_resis_Callback(hObject, eventdata, handles)
% hObject    handle to tor_resis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tor_resis as text
%        str2double(get(hObject,'String')) returns contents of tor_resis as a double


% --- Executes during object creation, after setting all properties.
function tor_resis_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tor_resis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function rpm_Callback(hObject, eventdata, handles)
% hObject    handle to rpm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of rpm as text
%        str2double(get(hObject,'String')) returns contents of rpm as a double


% --- Executes during object creation, after setting all properties.
function rpm_CreateFcn(hObject, eventdata, handles)
% hObject    handle to rpm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function numBlades_Callback(hObject, eventdata, handles)
% hObject    handle to numBlades (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of numBlades as text
%        str2double(get(hObject,'String')) returns contents of numBlades as a double


% --- Executes during object creation, after setting all properties.
function numBlades_CreateFcn(hObject, eventdata, handles)
% hObject    handle to numBlades (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function R_Callback(hObject, eventdata, handles)
% hObject    handle to R (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of R as text
%        str2double(get(hObject,'String')) returns contents of R as a double


% --- Executes during object creation, after setting all properties.
function R_CreateFcn(hObject, eventdata, handles)
% hObject    handle to R (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function hubR_Callback(hObject, eventdata, handles)
% hObject    handle to hubR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of hubR as text
%        str2double(get(hObject,'String')) returns contents of hubR as a double


% --- Executes during object creation, after setting all properties.
function hubR_CreateFcn(hObject, eventdata, handles)
% hObject    handle to hubR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function vo_design_Callback(hObject, eventdata, handles)
% hObject    handle to vo_design (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of vo_design as text
%        str2double(get(hObject,'String')) returns contents of vo_design as a double


% --- Executes during object creation, after setting all properties.
function vo_design_CreateFcn(hObject, eventdata, handles)
% hObject    handle to vo_design (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function pitch_Callback(hObject, eventdata, handles)
% hObject    handle to pitch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pitch as text
%        str2double(get(hObject,'String')) returns contents of pitch as a double


% --- Executes during object creation, after setting all properties.
function pitch_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pitch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function axes1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: place code in OpeningFcn to populate axes1


% --- Executes during object creation, after setting all properties.
function axes5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes5


% --- Executes during object creation, after setting all properties.
function axes6_CreateFcn(hObject2, eventdata, handles)
% hObject    handle to axes6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: place code in OpeningFcn to populate axes6



function rho_Callback(hObject, eventdata, handles)
% hObject    handle to rho (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of rho as text
%        str2double(get(hObject,'String')) returns contents of rho as a double


% --- Executes during object creation, after setting all properties.
function rho_CreateFcn(hObject, eventdata, handles)
% hObject    handle to rho (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function vis_Callback(hObject, eventdata, handles)
% hObject    handle to vis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of vis as text
%        str2double(get(hObject,'String')) returns contents of vis as a double


% --- Executes during object creation, after setting all properties.
function vis_CreateFcn(hObject, eventdata, handles)
% hObject    handle to vis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function E_Callback(hObject, eventdata, handles)
% hObject    handle to E (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of E as text
%        str2double(get(hObject,'String')) returns contents of E as a double


% --- Executes during object creation, after setting all properties.
function E_CreateFcn(hObject, eventdata, handles)
% hObject    handle to E (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function rh_Callback(hObject, eventdata, handles)
% hObject    handle to rh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of rh as text
%        str2double(get(hObject,'String')) returns contents of rh as a double


% --- Executes during object creation, after setting all properties.
function rh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to rh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ep_failure_Callback(hObject, eventdata, handles)
% hObject    handle to ep_failure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ep_failure as text
%        str2double(get(hObject,'String')) returns contents of ep_failure as a double


% --- Executes during object creation, after setting all properties.
function ep_failure_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ep_failure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function re_array_Callback(hObject, eventdata, handles)
% hObject    handle to re_array (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of re_array as text
%        str2double(get(hObject,'String')) returns contents of re_array as a double


% --- Executes during object creation, after setting all properties.
function re_array_CreateFcn(hObject, eventdata, handles)
% hObject    handle to re_array (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function MinsparTh_Callback(hObject, eventdata, handles)
% hObject    handle to MinsparTh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of MinsparTh as text
%        str2double(get(hObject,'String')) returns contents of MinsparTh as a double


% --- Executes during object creation, after setting all properties.
function MinsparTh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MinsparTh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function materialSF_Callback(hObject, eventdata, handles)
% hObject    handle to materialSF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of materialSF as text
%        str2double(get(hObject,'String')) returns contents of materialSF as a double


% --- Executes during object creation, after setting all properties.
function materialSF_CreateFcn(hObject, eventdata, handles)
% hObject    handle to materialSF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function loadSF_Callback(hObject, eventdata, handles)
% hObject    handle to loadSF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of loadSF as text
%        str2double(get(hObject,'String')) returns contents of loadSF as a double


% --- Executes during object creation, after setting all properties.
function loadSF_CreateFcn(hObject, eventdata, handles)
% hObject    handle to loadSF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function thick_Callback(hObject, eventdata, handles)
% hObject    handle to thick (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of thick as text
%        str2double(get(hObject,'String')) returns contents of thick as a double


% --- Executes during object creation, after setting all properties.
function thick_CreateFcn(hObject, eventdata, handles)
% hObject    handle to thick (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function L_Callback(hObject, eventdata, handles)
% hObject    handle to L (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of L as text
%        str2double(get(hObject,'String')) returns contents of L as a double


% --- Executes during object creation, after setting all properties.
function L_CreateFcn(hObject, eventdata, handles)
% hObject    handle to L (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function pop_Callback(hObject, eventdata, handles)
% hObject    handle to pop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pop as text
%        str2double(get(hObject,'String')) returns contents of pop as a double


% --- Executes during object creation, after setting all properties.
function pop_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function mingen_Callback(hObject, eventdata, handles)
% hObject    handle to mingen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mingen as text
%        str2double(get(hObject,'String')) returns contents of mingen as a double


% --- Executes during object creation, after setting all properties.
function mingen_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mingen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function val_Callback(hObject, eventdata, handles)
% hObject    handle to val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of val as text
%        str2double(get(hObject,'String')) returns contents of val as a double


% --- Executes during object creation, after setting all properties.
function val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
set(hObject, 'String', {'SBX','Uniform'});



function nm_Callback(hObject, eventdata, handles)
% hObject    handle to nm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nm as text
%        str2double(get(hObject,'String')) returns contents of nm as a double


% --- Executes during object creation, after setting all properties.
function nm_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function pm_Callback(hObject, eventdata, handles)
% hObject    handle to pm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pm as text
%        str2double(get(hObject,'String')) returns contents of pm as a double


% --- Executes during object creation, after setting all properties.
function pm_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function pc_Callback(hObject, eventdata, handles)
% hObject    handle to pc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pc as text
%        str2double(get(hObject,'String')) returns contents of pc as a double


% --- Executes during object creation, after setting all properties.
function pc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function nc_Callback(hObject, eventdata, handles)
% hObject    handle to nc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nc as text
%        str2double(get(hObject,'String')) returns contents of nc as a double


% --- Executes during object creation, after setting all properties.
function nc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function U_Callback(hObject, eventdata, handles)
% hObject    handle to U (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of U as text
%        str2double(get(hObject,'String')) returns contents of U as a double


% --- Executes during object creation, after setting all properties.
function U_CreateFcn(hObject, eventdata, handles)
% hObject    handle to U (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function uptwist_Callback(hObject, eventdata, handles)
% hObject    handle to uptwist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of uptwist as text
%        str2double(get(hObject,'String')) returns contents of uptwist as a double


% --- Executes during object creation, after setting all properties.
function uptwist_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uptwist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function dotwist_Callback(hObject, eventdata, handles)
% hObject    handle to dotwist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dotwist as text
%        str2double(get(hObject,'String')) returns contents of dotwist as a double


% --- Executes during object creation, after setting all properties.
function dotwist_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dotwist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function upchord_Callback(hObject, eventdata, handles)
% hObject    handle to upchord (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of upchord as text
%        str2double(get(hObject,'String')) returns contents of upchord as a double


% --- Executes during object creation, after setting all properties.
function upchord_CreateFcn(hObject, eventdata, handles)
% hObject    handle to upchord (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function dochord_Callback(hObject, eventdata, handles)
% hObject    handle to dochord (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dochord as text
%        str2double(get(hObject,'String')) returns contents of dochord as a double


% --- Executes during object creation, after setting all properties.
function dochord_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dochord (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function N_Callback(hObject, eventdata, handles)
% hObject    handle to N (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of N as text
%        str2double(get(hObject,'String')) returns contents of N as a double


% --- Executes during object creation, after setting all properties.
function N_CreateFcn(hObject, eventdata, handles)
% hObject    handle to N (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function airfoildata_Callback(hObject, eventdata, handles)
% hObject    handle to airfoildata (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of airfoildata as text
%        str2double(get(hObject,'String')) returns contents of airfoildata as a double


% --- Executes during object creation, after setting all properties.
function airfoildata_CreateFcn(hObject, eventdata, handles)
% hObject    handle to airfoildata (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function airfoilprofile_Callback(hObject, eventdata, handles)
% hObject    handle to airfoilprofile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of airfoilprofile as text
%        str2double(get(hObject,'String')) returns contents of airfoilprofile as a double


% --- Executes during object creation, after setting all properties.
function airfoilprofile_CreateFcn(hObject, eventdata, handles)
% hObject    handle to airfoilprofile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Ve50_Callback(hObject, eventdata, handles)
% hObject    handle to Ve50 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Ve50 as text
%        str2double(get(hObject,'String')) returns contents of Ve50 as a double


% --- Executes during object creation, after setting all properties.
function Ve50_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Ve50 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function pvalstruct_Callback(hObject, eventdata, handles)
% hObject    handle to pvalstruct (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pvalstruct as text
%        str2double(get(hObject,'String')) returns contents of pvalstruct as a double


% --- Executes during object creation, after setting all properties.
function pvalstruct_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pvalstruct (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function pval_Callback(hObject, eventdata, handles)
% hObject    handle to pval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pval as text
%        str2double(get(hObject,'String')) returns contents of pval as a double


% --- Executes during object creation, after setting all properties.
function pval_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function valBplot_Callback(hObject, eventdata, handles)
% hObject    handle to valBplot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of valBplot as text
%        str2double(get(hObject,'String')) returns contents of valBplot as a double


% --- Executes during object creation, after setting all properties.
function valBplot_CreateFcn(hObject, eventdata, handles)
% hObject    handle to valBplot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function z0_Callback(hObject, eventdata, handles)
% hObject    handle to z0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of z0 as text
%        str2double(get(hObject,'String')) returns contents of z0 as a double


% --- Executes during object creation, after setting all properties.
function z0_CreateFcn(hObject, eventdata, handles)
% hObject    handle to z0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Ratio_Callback(hObject, eventdata, handles)
% hObject    handle to Ratio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Ratio as text
%        str2double(get(hObject,'String')) returns contents of Ratio as a double


% --- Executes during object creation, after setting all properties.
function Ratio_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Ratio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Htower_Callback(hObject, eventdata, handles)
% hObject    handle to Htower (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Htower as text
%        str2double(get(hObject,'String')) returns contents of Htower as a double


% --- Executes during object creation, after setting all properties.
function Htower_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Htower (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function a_top_Callback(hObject, eventdata, handles)
% hObject    handle to a_top (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of a_top as text
%        str2double(get(hObject,'String')) returns contents of a_top as a double


% --- Executes during object creation, after setting all properties.
function a_top_CreateFcn(hObject, eventdata, handles)
% hObject    handle to a_top (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function a_ground_Callback(hObject, eventdata, handles)
% hObject    handle to a_ground (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of a_ground as text
%        str2double(get(hObject,'String')) returns contents of a_ground as a double


% --- Executes during object creation, after setting all properties.
function a_ground_CreateFcn(hObject, eventdata, handles)
% hObject    handle to a_ground (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Gamma_Callback(hObject, eventdata, handles)
% hObject    handle to Gamma (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Gamma as text
%        str2double(get(hObject,'String')) returns contents of Gamma as a double


% --- Executes during object creation, after setting all properties.
function Gamma_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Gamma (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function shaftlength_Callback(hObject, eventdata, handles)
% hObject    handle to shaftlength (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of shaftlength as text
%        str2double(get(hObject,'String')) returns contents of shaftlength as a double


% --- Executes during object creation, after setting all properties.
function shaftlength_CreateFcn(hObject, eventdata, handles)
% hObject    handle to shaftlength (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Tscale_Callback(hObject, eventdata, handles)
% hObject    handle to Tscale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Tscale as text
%        str2double(get(hObject,'String')) returns contents of Tscale as a double


% --- Executes during object creation, after setting all properties.
function Tscale_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Tscale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Tinten_Callback(hObject, eventdata, handles)
% hObject    handle to Tinten (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Tinten as text
%        str2double(get(hObject,'String')) returns contents of Tinten as a double


% --- Executes during object creation, after setting all properties.
function Tinten_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Tinten (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function hblunt_Callback(hObject, eventdata, handles)
% hObject    handle to hblunt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of hblunt as text
%        str2double(get(hObject,'String')) returns contents of hblunt as a double


% --- Executes during object creation, after setting all properties.
function hblunt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to hblunt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function r0_Callback(hObject, eventdata, handles)
% hObject    handle to r0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of r0 as text
%        str2double(get(hObject,'String')) returns contents of r0 as a double


% --- Executes during object creation, after setting all properties.
function r0_CreateFcn(hObject, eventdata, handles)
% hObject    handle to r0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function round_Callback(hObject, eventdata, handles)
% hObject    handle to round (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of round as text
%        str2double(get(hObject,'String')) returns contents of round as a double


% --- Executes during object creation, after setting all properties.
function round_CreateFcn(hObject, eventdata, handles)
% hObject    handle to round (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

set(hObject, 'String', {'Squared','Rounded'});



function h0_Callback(hObject, eventdata, handles)
% hObject    handle to h0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of h0 as text
%        str2double(get(hObject,'String')) returns contents of h0 as a double


% --- Executes during object creation, after setting all properties.
function h0_CreateFcn(hObject, eventdata, handles)
% hObject    handle to h0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function TEangle_Callback(hObject, eventdata, handles)
% hObject    handle to TEangle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of TEangle as text
%        str2double(get(hObject,'String')) returns contents of TEangle as a double


% --- Executes during object creation, after setting all properties.
function TEangle_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TEangle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function PSI_Callback(hObject, eventdata, handles)
% hObject    handle to PSI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of PSI as text
%        str2double(get(hObject,'String')) returns contents of PSI as a double


% --- Executes during object creation, after setting all properties.
function PSI_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PSI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function add_noise_Callback(hObject, eventdata, handles)
% hObject    handle to add_noise (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of add_noise
if (get(hObject,'Value') == get(hObject,'Max'))
    % Checkbox is checked-take approriate action
    set(handles.zo,'Enable','on');
    set(handles.a_top,'Enable','on');
    set(handles.a_ground,'Enable','on');
    set(handles.Htower,'Enable','on');
    set(handles.shaftlength,'Enable','on');
    set(handles.Gamma,'Enable','on');
    set(handles.estimate,'Enable','on');
    set(handles.Tscale,'Enable','on');
    set(handles.Tinten,'Enable','on');
    set(handles.round,'Enable','on');
    set(handles.hblunt,'Enable','on');
    set(handles.Ratio,'Enable','on');
    set(handles.TEangle,'Enable','on');
    set(handles.h0,'Enable','on');
    set(handles.PSI,'Enable','on');
    set(handles.noisedatabase,'Enable','on');
    set(handles.r0,'Enable','on');
    set(handles.w4,'Enable','on');
    set(handles.Rref,'Enable','on');
    set(handles.Aref,'Enable','on');
else
    % Checkbox is not checked-take approriate action
    set(handles.zo,'Enable','off');
    set(handles.a_top,'Enable','off');
    set(handles.a_ground,'Enable','off');
    set(handles.Htower,'Enable','off');
    set(handles.shaftlength,'Enable','off');
    set(handles.Gamma,'Enable','off');
    set(handles.estimate,'Enable','off');
    set(handles.Tscale,'Enable','off');
    set(handles.Tinten,'Enable','off');
    set(handles.round,'Enable','off');
    set(handles.hblunt,'Enable','off');
    set(handles.Ratio,'Enable','off');
    set(handles.TEangle,'Enable','off');
    set(handles.h0,'Enable','off');
    set(handles.PSI,'Enable','off');
    set(handles.noisedatabase,'Enable','off');
    set(handles.r0,'Enable','off');
    set(handles.w4,'Enable','off');
    set(handles.w4, 'String', num2str(0));
    set(handles.Rref,'Enable','off');
    set(handles.Aref,'Enable','off');
end

guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function add_noise_CreateFcn(hObject, eventdata, handles)
% hObject    handle to add_noise (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function w1_Callback(hObject, eventdata, handles)
% hObject    handle to w1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of w1 as text
%        str2double(get(hObject,'String')) returns contents of w1 as a double


% --- Executes during object creation, after setting all properties.
function w1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to w1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function w2_Callback(hObject, eventdata, handles)
% hObject    handle to w2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of w2 as text
%        str2double(get(hObject,'String')) returns contents of w2 as a double


% --- Executes during object creation, after setting all properties.
function w2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to w2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function w3_Callback(hObject, eventdata, handles)
% hObject    handle to w3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of w3 as text
%        str2double(get(hObject,'String')) returns contents of w3 as a double


% --- Executes during object creation, after setting all properties.
function w3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to w3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function w4_Callback(hObject, eventdata, handles)
% hObject    handle to w4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of w4 as text
%        str2double(get(hObject,'String')) returns contents of w4 as a double


% --- Executes during object creation, after setting all properties.
function w4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to w4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in uipanel124.
function loadcases_Callback(hObject, eventdata, handles)
% hObject    handle to uipanel124 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns uipanel124 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from uipanel124
thevalue=get(hObject,'Value');
if (thevalue == 1)||(thevalue == 2)
    % Checkbox is checked-take approriate action
    set(handles.rotoroverspeed,'Enable','off');
else
    % Checkbox is not checked-take approriate action
    set(handles.rotoroverspeed,'Enable','on');
end

% --- Executes during object creation, after setting all properties.
function loadcases_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uipanel124 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

set(hObject, 'String', {'---Load Cases---','Parked Condition (default)','Rotor Overspeed'});



function rotoroverspeed_Callback(hObject, eventdata, handles)
% hObject    handle to rotoroverspeed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of rotoroverspeed as text
%        str2double(get(hObject,'String')) returns contents of rotoroverspeed as a double


% --- Executes during object creation, after setting all properties.
function rotoroverspeed_CreateFcn(hObject, eventdata, handles)
% hObject    handle to rotoroverspeed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
set(hObject,'Enable','off');


% --- Executes on button press in runBaseline.
function runBaseline_Callback(hObject, eventdata, handles)
% hObject    handle to runBaseline (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of runBaseline
if get(hObject,'Value')==1
    set(handles.turbdata,'Enable','on')
elseif get(hObject,'Value')==0
    set(handles.turbdata,'Enable','off')
end

% --- Executes during object creation, after setting all properties.
function runBaseline_CreateFcn(hObject, eventdata, handles)
% hObject    handle to add_noise (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
if get(hObject,'Value')==1
    set(handles.turbdata,'Enable','on')
elseif get(hObject,'Value')==0
    set(handles.turbdata,'Enable','off')
end

% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%Airfoil profile, lift and drag data
try
    data.airfoildatastring=(get(handles.airfoildata,'String'));
    data.profiledatastring=(get(handles.airfoilprofile,'String'));
    data.re_array=str2num(get(handles.re_array,'String'));
    [h1] = plotairfoildata(data);
catch
    disp('----!!!!!PLEASE CHECK YOUR AIRFOIL DATA!!!!----')
    return
end


% --- Executes on selection change in spacing.
function spacing_Callback(hObject, eventdata, handles)
% hObject    handle to spacing (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns spacing contents as cell array
%        contents{get(hObject,'Value')} returns selected item from spacing


% --- Executes during object creation, after setting all properties.
function spacing_CreateFcn(hObject, eventdata, handles)
% hObject    handle to spacing (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

set(hObject, 'String', {'---Control Point Spacing---','uniform (default)','half-cosine (finer towards hub)','half-cosine (finer towards tip)','full-cosine'});



function turbdata_Callback(hObject, eventdata, handles)
% hObject    handle to turbdata (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of turbdata as text
%        str2double(get(hObject,'String')) returns contents of turbdata as a double


% --- Executes during object creation, after setting all properties.
function turbdata_CreateFcn(hObject, eventdata, handles)
% hObject    handle to turbdata (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in optprogress.
function optprogress_Callback(hObject, eventdata, handles)
% hObject    handle to optprogress (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of optprogress


% --- Executes on button press in flag3.
function flag3_Callback(hObject, eventdata, handles)
% hObject    handle to flag3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of flag3


function relax_Callback(hObject, eventdata, handles)
% hObject    handle to relax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of relax as text
%        str2double(get(hObject,'String')) returns contents of relax as a double


% --- Executes during object creation, after setting all properties.
function relax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to relax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in printcoords.
function printcoords_Callback(hObject, eventdata, handles)
% hObject    handle to printcoords (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of printcoords


% --- Executes on selection change in spacingBE.
function spacingBE_Callback(hObject, eventdata, handles)
% hObject    handle to spacingBE (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns spacingBE contents as cell array
%        contents{get(hObject,'Value')} returns selected item from spacingBE


% --- Executes during object creation, after setting all properties.
function spacingBE_CreateFcn(hObject, eventdata, handles)
% hObject    handle to spacingBE (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

set(hObject, 'String', {'---Blade Element Spacing---','half-cosine (finer towards tip) (default)','half-cosine (finer towards hub)','full-cosine','uniform'});


% --- Executes on selection change in zo.
function zo_Callback(hObject, eventdata, handles)
% hObject    handle to zo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns zo contents as cell array
%        contents{get(hObject,'Value')} returns selected item from zo


% --- Executes during object creation, after setting all properties.
function zo_CreateFcn(hObject, eventdata, handles)
% hObject    handle to zo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

set(hObject, 'String', {'Very smooth, ice or mud (0.01)','Calm open sea (0.20)','Blown sea (0.50)','Snow surface (3.00)',...
    'Lawn grass (8.00)','Rough pasture (10.00)','Fallow field (30.00)','Crops (50.00)','Few trees (100.00)',...
    'Many trees, hedges (250.00)','Forest and woodlands (500.00)','Suburbs (1500.00)',...
    'Centers of cities with tall buildings (3000.00)'});


% --- Executes on button press in observerpos.
function observerpos_Callback(hObject, eventdata, handles)
% hObject    handle to observerpos (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Display image file
k1=figure(50);
observerpos = imread('imagesGUI\observerpos.jpg');
imagesc(observerpos);
set(gca,'visible','off')
set(gca,'LooseInset',get(gca,'TightInset'))
set(k1, 'units', 'centimeters', 'pos', [2.0 2.0 18.1184 14.1112])


% --- Executes on button press in estimate.
function estimate_Callback(hObject, eventdata, handles)
% hObject    handle to estimate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of estimate


% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%Display image file
k1=figure(50);
gamma = imread('imagesGUI\gamma.jpg');
imagesc(gamma);
set(gca,'visible','off')
set(gca,'LooseInset',get(gca,'TightInset'))
set(k1, 'units', 'centimeters', 'pos', [2.0 2.0 29.422 11.994])


% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%Display image file
k1=figure(50);
blunt = imread('imagesGUI\blunt.jpg');
imagesc(blunt);
set(gca,'visible','off')
set(gca,'LooseInset',get(gca,'TightInset'))
set(k1, 'units', 'centimeters', 'pos', [2.0 2.0 18.732 5.715])


% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%Display image file
k1=figure(50);
tipfactor = imread('imagesGUI\tipfactor.jpg');
imagesc(tipfactor);
set(gca,'visible','off')
set(gca,'LooseInset',get(gca,'TightInset'))
set(k1, 'units', 'centimeters', 'pos', [2.0 2.0 28.998 16.686])


% --- Executes on button press in xfoilexecute.
function xfoilexecute_Callback(hObject, eventdata, handles)
% hObject    handle to xfoilexecute (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of xfoilexecute


% --- Executes on selection change in Ncrit.
function Ncrit_Callback(hObject, eventdata, handles)
% hObject    handle to Ncrit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns Ncrit contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Ncrit
NcritVal=get(handles.Ncrit,'Value');
if NcritVal==3
    set(handles.trip, 'String', 'tripped');
else
    set(handles.trip, 'String', 'untripped');
end

% --- Executes during object creation, after setting all properties.
function Ncrit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Ncrit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

set(hObject, 'String', {'---Transition trigger (N_crit)---','untriggered: N_crit = 9 (default)','triggered: N_crit = 4'});



function ITER_Callback(hObject, eventdata, handles)
% hObject    handle to ITER (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ITER as text
%        str2double(get(hObject,'String')) returns contents of ITER as a double


% --- Executes during object creation, after setting all properties.
function ITER_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ITER (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton8.
function pushbutton8_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%Display image file
k1=figure(50);
tipfactor = imread('imagesGUI\ncrit.jpg');
imagesc(tipfactor);
set(gca,'visible','off')
set(gca,'LooseInset',get(gca,'TightInset'))
set(k1, 'units', 'centimeters', 'pos', [2.0 2.0 26.603 16.686])


% --- Executes on button press in generateNoise.
function generateNoise_Callback(hObject, eventdata, handles)
% hObject    handle to generateNoise (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%XFOIL Data

%XFOIL data
itersXFOIL=str2num(get(handles.ITER,'String'));
NcritVal=get(handles.Ncrit,'Value');
datatemp.airfoildatastring=(get(handles.airfoildata,'String'));
datatemp.profiledatastring=(get(handles.airfoilprofile,'String'));
[~] = airfoildata(datatemp);

%Retrieve filename and path from user
[noise_FileName,noise_PathName] = uiputfile([pwd '\noisedatabase\*.dat'],'Save generated noise database');
if ~ischar(noise_FileName)
    return; %if user canceled exit this callback
end
%construct the path name of the save location
saveNoiseDataName = fullfile(noise_PathName,noise_FileName); 

Ncrit=9;
%logi_trip --> Condition for boundary layer tripping: tripped=0, untripped=1, partially tripped=2
if NcritVal==3
    Ncrit=4;
end
Rref=str2num(get(handles.Rref,'String'));  % Reference Reynolds number
Aref=str2num(get(handles.Aref,'String'));  % Reference alfa number

[~]=xfoil_execution(Rref,Aref,itersXFOIL,Ncrit,saveNoiseDataName);


% --- Executes on button press in Load_GUI.
function Load_GUI_Callback(hObject, eventdata, handles)
% hObject    handle to Load_GUI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[LoadGUI_FileName, LoadGUI_PathName] = uigetfile([pwd '\GUI_settings\*.fig'],'Load GUI settings');

if ~ischar(LoadGUI_FileName)
    return; %if user canceled exit this callback
end

%closes the old gui
close(gcf);

%construct the path name of the file to be loaded
loadDataName = fullfile(LoadGUI_PathName,LoadGUI_FileName);

%load the settings, which creates a new gui
f2=hgload(loadDataName);
movegui(f2,'center') 

% --- Executes on button press in Save_GUI.
function Save_GUI_Callback(hObject, eventdata, handles)
% hObject    handle to Save_GUI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[GUI_FileName,GUI_PathName] = uiputfile([pwd '\GUI_settings\*.fig'],'Save GUI settings');

if ~ischar(GUI_FileName)
    return; %if user canceled exit this callback
end

%construct the path name of the save location
saveGUIDataName = fullfile(GUI_PathName,GUI_FileName); 
 
%saves the gui data
hgsave(saveGUIDataName);


% --- Executes on button press in pushbutton12.
function pushbutton12_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
LoadAirfoilData_FileName = uigetfile([pwd '\airfoildatabase\*.txt'],'Load airfoil data - FILE MUST BE LOCATED IN ''airfoildatabase'' FOLDER');

if ~ischar(LoadAirfoilData_FileName)
    return; %if user canceled exit this callback
end

%construct the path name of the file to be loaded
loadAirfoilDataName = fullfile('airfoildatabase\',LoadAirfoilData_FileName);

set(handles.airfoildata,'String',loadAirfoilDataName);




% --- Executes on button press in pushbutton14.
function pushbutton14_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
LoadAirfoilProfile_FileName = uigetfile([pwd '\airfoildatabase\*.txt'],'Load airfoil profile - FILE MUST BE LOCATED IN ''airfoildatabase'' FOLDER');

if ~ischar(LoadAirfoilProfile_FileName)
    return; %if user canceled exit this callback
end

%construct the path name of the file to be loaded
loadAirfoilProfileName = fullfile('airfoildatabase\',LoadAirfoilProfile_FileName);

set(handles.airfoilprofile,'String',loadAirfoilProfileName);


% --- Executes on button press in pushbutton15.
function pushbutton15_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[LoadTurbData_FileName LoadTurbData_PathName] = uigetfile([pwd '\turbinedatabase\*.txt'],'Load turbine data');

if ~ischar(LoadTurbData_FileName)
    return; %if user canceled exit this callback
end

%construct the path name of the file to be loaded
loadTurbDataName = fullfile(LoadTurbData_PathName,LoadTurbData_FileName);

set(handles.turbdata,'String',loadTurbDataName);


% --- Executes on button press in test_baseline.
function test_baseline_Callback(hObject, eventdata, handles)
% hObject    handle to test_baseline (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if (get(handles.runBaseline,'Value')==1)
    [data] = collect_data(handles);
end


% --- Executes on button press in pushbutton17.
function pushbutton17_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[LoadReArray_FileName LoadReArray_PathName] = uigetfile([pwd '\airfoildatabase\*.txt'],'Load Reynolds number array');

if ~ischar(LoadReArray_FileName)
    return; %if user canceled exit this callback
end

%construct the path name of the file to be loaded
loadReArrayName = fullfile(LoadReArray_PathName,LoadReArray_FileName);

array_re=importdata(loadReArrayName);

set(handles.re_array,'String',['[' num2str(array_re) ']']);


% --- Executes on button press in parallel.
function parallel_Callback(hObject, eventdata, handles)
% hObject    handle to parallel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of parallel


% --- Executes on button press in cptsr.
function cptsr_Callback(hObject, eventdata, handles)
% hObject    handle to cptsr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cptsr



function tsrmax_Callback(hObject, eventdata, handles)
% hObject    handle to tsrmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tsrmax as text
%        str2double(get(hObject,'String')) returns contents of tsrmax as a double


% --- Executes during object creation, after setting all properties.
function tsrmax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tsrmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function tsrmin_Callback(hObject, eventdata, handles)
% hObject    handle to tsrmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tsrmin as text
%        str2double(get(hObject,'String')) returns contents of tsrmin as a double


% --- Executes during object creation, after setting all properties.
function tsrmin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tsrmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton18.
function pushbutton18_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton18 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
h120 = msgbox('Cp - TSR curve will only be exported into Excel for the baseline and final blade designs.','NOTE','help');



function noisedatabase_Callback(hObject, eventdata, handles)
% hObject    handle to noisedatabase (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of noisedatabase as text
%        str2double(get(hObject,'String')) returns contents of noisedatabase as a double


% --- Executes during object creation, after setting all properties.
function noisedatabase_CreateFcn(hObject, eventdata, handles)
% hObject    handle to noisedatabase (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton19.
function pushbutton19_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[LoadNoiseData_FileName LoadNoiseData_PathName] = uigetfile([pwd '\noisedatabase\*.dat'],'Load airfoil noise database');

if ~ischar(LoadNoiseData_FileName)
    return; %if user canceled exit this callback
end

%construct the path name of the file to be loaded
loadNoiseDataName = fullfile(LoadNoiseData_PathName,LoadNoiseData_FileName);

set(handles.noisedatabase,'String',loadNoiseDataName);



function lambda_design_Callback(hObject, eventdata, handles)
% hObject    handle to lambda_design (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of lambda_design as text
%        str2double(get(hObject,'String')) returns contents of lambda_design as a double


% --- Executes during object creation, after setting all properties.
function lambda_design_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lambda_design (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in attention.
function attention_Callback(hObject, eventdata, handles)
% hObject    handle to attention (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
h121 = msgbox('Prior to running the optimization, you must set up the number of cores using the "matlabpool open #" command.','NOTE','help');


% --- Executes on button press in attention2.
function attention2_Callback(hObject, eventdata, handles)
% hObject    handle to attention2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
h122 = msgbox('Note: w1 + w2 + ... = 1','NOTE','help');


% --- Executes on button press in optAEP.
function optAEP_Callback(hObject, eventdata, handles)
% hObject    handle to optAEP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of optAEP
if (get(hObject,'Value') == get(hObject,'Max'))
    set(handles.probtype,'Value',1)
    set(handles.probtype,'Enable','on')
    set(handles.k_weibull,'Enable','off');
    set(handles.A_weibull,'Enable','off');
    set(handles.vcutin,'Enable','on');
    set(handles.vcutout,'Enable','on');
    set(handles.mws,'Enable','on');
    set(handles.winddata,'Enable','off');
        
    set(handles.cptext, 'String', 'AEP');
    set(handles.maxpower,'Enable','on');
else
    % Checkbox is not checked-take approriate action
    set(handles.probtype,'Value',1)
    set(handles.probtype,'Enable','off')
    set(handles.k_weibull,'Enable','off');
    set(handles.A_weibull,'Enable','off');
    set(handles.vcutin,'Enable','off');
    set(handles.vcutout,'Enable','off');
    set(handles.mws,'Enable','off');
    set(handles.winddata,'Enable','off');

    set(handles.cptext, 'String', 'Cp');
    set(handles.maxpower,'Enable','off');
end


function winddata_Callback(hObject, eventdata, handles)
% hObject    handle to winddata (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of winddata as text
%        str2double(get(hObject,'String')) returns contents of winddata as a double


% --- Executes during object creation, after setting all properties.
function winddata_CreateFcn(hObject, eventdata, handles)
% hObject    handle to winddata (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
set(hObject,'Enable','off');


% --- Executes on button press in pushbutton22.
function pushbutton22_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[LoadWindData_FileName LoadWindData_PathName] = uigetfile([pwd '\wind_data\*.dat'],'Load wind data');

if ~ischar(LoadWindData_FileName)
    return; %if user canceled exit this callback
end

%construct the path name of the file to be loaded
loadWindDataName = fullfile(LoadWindData_PathName,LoadWindData_FileName);

set(handles.winddata,'String',loadWindDataName);

% --- Executes on button press in add_starting.
function add_starting_Callback(hObject, eventdata, handles)
% hObject    handle to add_starting (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of add_starting
if (get(hObject,'Value') == get(hObject,'Max'))
    % Checkbox is checked-take approriate action
    set(handles.u_start,'Enable','on');
    set(handles.lambda_start,'Enable','on');
    set(handles.Jgen,'Enable','on');
    set(handles.tor_resis,'Enable','on');
    set(handles.w2,'Enable','on');
else
    % Checkbox is not checked-take approriate action
    set(handles.u_start,'Enable','off');
    set(handles.lambda_start,'Enable','off');
    set(handles.Jgen,'Enable','off');
    set(handles.tor_resis,'Enable','off');
    set(handles.w2,'Enable','off');
    set(handles.w2, 'String', num2str(0));
end


% --- Executes during object creation, after setting all properties.
function add_structure_CreateFcn(hObject, eventdata, handles)
% hObject    handle to add_structure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
set(hObject,'String',{'Solid','Hollow'});

% --- Executes on button press in add_structure.
function add_structure_Callback(hObject, eventdata, handles)
% hObject    handle to add_starting (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of add_starting
if (get(hObject,'Value') == 2)
    % Checkbox is checked-take approriate action
    set(handles.MinsparTh,'Enable','on');
else
    % Checkbox is not checked-take approriate action
    set(handles.MinsparTh,'Enable','off');
end

% --- Executes on button press in perform_structure.
function perform_structure_Callback(hObject, eventdata, handles)
% hObject    handle to perform_structure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of perform_structure
if (get(hObject,'Value') == 1)
    % Checkbox is checked-take approriate action
    set(handles.E,'Enable','on');
    set(handles.materialSF,'Enable','on');
    set(handles.rh,'Enable','on');
    set(handles.U,'Enable','on');
    set(handles.L,'Enable','on');
    set(handles.ep_failure,'Enable','on');
    if handles.add_structure==2
        set(handles.MinsparTh,'Enable','on');
    end
    set(handles.w3,'Enable','on');
    set(handles.add_structure,'Enable','on');
else
    % Checkbox is not checked-take approriate action
    set(handles.E,'Enable','off');
    set(handles.materialSF,'Enable','off');
    set(handles.rh,'Enable','off');
    set(handles.U,'Enable','off');
    set(handles.L,'Enable','off');
    set(handles.ep_failure,'Enable','off');
    set(handles.MinsparTh,'Enable','off');
    set(handles.w3, 'String', num2str(0));
    set(handles.w3,'Enable','off');
    set(handles.add_structure, 'Value', 1);
    set(handles.add_structure,'Enable','off');
end


% --- Executes on button press in check_noise.
function check_noise_Callback(hObject, eventdata, handles)
% hObject    handle to check_noise (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
noisedatabase=(get(handles.noisedatabase,'String'));
try
    noisedata_dat=importdata(noisedatabase);
catch
    disp('Please input a correct airfoil noise database file.')
    return
end
no_RE=size(noisedata_dat,1)/2;
P=noisedata_dat(1:no_RE,:);
S=noisedata_dat((no_RE+1):end,:);

Rref=str2num(get(handles.Rref,'String'));
Aref=str2num(get(handles.Aref,'String'));

if size(P,1)~=length(Rref) || size(P,2)~=length(Aref)
    disp('Input noise database file does not match reference alfa and/or Reynolds number range(s)')
    return
end

h6=figure(50);
subplot(1,2,1);
[X Y]=meshgrid(Aref,Rref/10^6);
surf(X,Y,P)
xlabel('AoA (degrees)')
ylabel('Re (10^6)')
zlabel('\delta^* at TE (m)')
title('\delta^* at TE - Pressure Side')

subplot(1,2,2);
surf(X,Y,S)
xlabel('AoA (degrees)')
ylabel('Re (10^6)')
zlabel('\delta^* at TE (m)')
title('\delta^* at TE - Suction Side')

set(h6, 'units', 'centimeters', 'pos', [0 2 25 12])


function maxpower_Callback(hObject, eventdata, handles)
% hObject    handle to maxpower (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of maxpower as text
%        str2double(get(hObject,'String')) returns contents of maxpower as a double


% --- Executes during object creation, after setting all properties.
function maxpower_CreateFcn(hObject, eventdata, handles)
% hObject    handle to maxpower (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in popsizemusteven.
function popsizemusteven_Callback(hObject, eventdata, handles)
% hObject    handle to popsizemusteven (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 h123 = msgbox('Population size must be even when divided by 2','NOTE','help');

 

% --- Executes on button press in SaveResults.
function SaveResults_Callback(hObject, eventdata, handles)
% hObject    handle to SaveResults (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[Results_FileName,Results_PathName] = uiputfile([pwd '\output\*.xls'],'Save results file');

if ~ischar(Results_FileName)
    return; %if user canceled exit this callback
end

%construct the path name of the save location
saveResultsDataName = fullfile(Results_PathName,Results_FileName); 
 
%Store string
set(handles.saveTextbox,'String',saveResultsDataName);



function saveTextbox_Callback(hObject, eventdata, handles)
% hObject    handle to saveTextbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of saveTextbox as text
%        str2double(get(hObject,'String')) returns contents of saveTextbox as a double


% --- Executes during object creation, after setting all properties.
function saveTextbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to saveTextbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in probtype.
function probtype_Callback(hObject, eventdata, handles)
% hObject    handle to probtype (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns probtype contents as cell array
%        contents{get(hObject,'Value')} returns selected item from probtype
thevalue=get(hObject,'Value');
if thevalue == 3
    set(handles.k_weibull,'Enable','on');
    set(handles.A_weibull,'Enable','on');
    set(handles.vcutin,'Enable','on');
    set(handles.vcutout,'Enable','on');
    set(handles.mws,'Enable','off');
    set(handles.winddata,'Enable','off');
elseif thevalue == 4
    set(handles.k_weibull,'Enable','off');
    set(handles.A_weibull,'Enable','off');
    set(handles.vcutin,'Enable','on');
    set(handles.vcutout,'Enable','on');
    set(handles.mws,'Enable','off');
    set(handles.winddata,'Enable','on');
else
    set(handles.k_weibull,'Enable','off');
    set(handles.A_weibull,'Enable','off');
    set(handles.vcutin,'Enable','on');
    set(handles.vcutout,'Enable','on');
    set(handles.mws,'Enable','on');
    set(handles.winddata,'Enable','off');
end

% --- Executes during object creation, after setting all properties.
function probtype_CreateFcn(hObject, eventdata, handles)
% hObject    handle to probtype (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
set(hObject,'String',{'---Probability type---','Rayleigh (default)','Weibull','User-defined'});
set(hObject,'Enable','off');

function vcutin_Callback(hObject, eventdata, handles)
% hObject    handle to vcutin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of vcutin as text
%        str2double(get(hObject,'String')) returns contents of vcutin as a double


% --- Executes during object creation, after setting all properties.
function vcutin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to vcutin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
set(hObject,'Enable','off');


function vcutout_Callback(hObject, eventdata, handles)
% hObject    handle to vcutout (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of vcutout as text
%        str2double(get(hObject,'String')) returns contents of vcutout as a double


% --- Executes during object creation, after setting all properties.
function vcutout_CreateFcn(hObject, eventdata, handles)
% hObject    handle to vcutout (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
set(hObject,'Enable','off');


function mws_Callback(hObject, eventdata, handles)
% hObject    handle to mws (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mws as text
%        str2double(get(hObject,'String')) returns contents of mws as a double


% --- Executes during object creation, after setting all properties.
function mws_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mws (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
set(hObject,'Enable','off');


function A_weibull_Callback(hObject, eventdata, handles)
% hObject    handle to A_weibull (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of A_weibull as text
%        str2double(get(hObject,'String')) returns contents of A_weibull as a double


% --- Executes during object creation, after setting all properties.
function A_weibull_CreateFcn(hObject, eventdata, handles)
% hObject    handle to A_weibull (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
set(hObject,'Enable','off');


function k_weibull_Callback(hObject, eventdata, handles)
% hObject    handle to k_weibull (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of k_weibull as text
%        str2double(get(hObject,'String')) returns contents of k_weibull as a double


% --- Executes during object creation, after setting all properties.
function k_weibull_CreateFcn(hObject, eventdata, handles)
% hObject    handle to k_weibull (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
set(hObject,'Enable','off');


% --- Executes on button press in pushbutton27.
function pushbutton27_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton27 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
            %User-Defined Probability Distribution
 
%AEP data
tempdata.optAEP=get(handles.optAEP,'Value');
tempdata.winddata=get(handles.winddata,'String');
tempdata.mws=str2num(get(handles.mws,'String'));
tempdata.probtype=get(handles.probtype,'Value');
tempdata.vcutin=str2num(get(handles.vcutin,'String'));
tempdata.vcutout=str2num(get(handles.vcutout,'String'));
tempdata.k_weibull=str2num(get(handles.k_weibull,'String'));
tempdata.A_weibull=str2num(get(handles.A_weibull,'String'));

[AEP]=aep(0,tempdata,1);



function Aref_Callback(hObject, eventdata, handles)
% hObject    handle to Aref (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Aref as text
%        str2double(get(hObject,'String')) returns contents of Aref as a double


% --- Executes during object creation, after setting all properties.
function Aref_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Aref (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton29.
function pushbutton29_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton29 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[LoadRref_FileName LoadRref_PathName] = uigetfile([pwd '\noisedatabase\*.dat'],'Load reference Reynolds numbers for noise database');

if ~ischar(LoadRref_FileName)
    return; %if user canceled exit this callback
end

%construct the path name of the file to be loaded
loadRrefName = fullfile(LoadRref_PathName,LoadRref_FileName);

Rref=importdata(loadRrefName);

set(handles.Rref,'String',['[' num2str(Rref.data') ']']);


function Rref_Callback(hObject, eventdata, handles)
% hObject    handle to Rref (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Rref as text
%        str2double(get(hObject,'String')) returns contents of Rref as a double


% --- Executes during object creation, after setting all properties.
function Rref_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Rref (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton30.
function pushbutton30_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton30 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[LoadAref_FileName LoadAref_PathName] = uigetfile([pwd '\noisedatabase\*.dat'],'Load reference angle of attack for noise database');

if ~ischar(LoadAref_FileName)
    return; %if user canceled exit this callback
end

%construct the path name of the file to be loaded
loadArefName = fullfile(LoadAref_PathName,LoadAref_FileName);

Aref=importdata(loadArefName);

set(handles.Aref,'String',['[' num2str(Aref.data') ']']);



function numcores_Callback(hObject, eventdata, handles)
% hObject    handle to numcores (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of numcores as text
%        str2double(get(hObject,'String')) returns contents of numcores as a double


% --- Executes during object creation, after setting all properties.
function numcores_CreateFcn(hObject, eventdata, handles)
% hObject    handle to numcores (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
