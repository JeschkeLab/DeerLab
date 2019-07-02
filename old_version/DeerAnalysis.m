function varargout = DeerAnalysis(varargin)
% DEERANALYSIS M-file for DeerAnalysis.fig
%      DEERANALYSIS, by itself, creates a new DEERANALYSIS or raises the existing
%      singleton*.
%
%      H = DEERANALYSIS returns the handle to a new DEERANALYSIS or the handle to
%      the existing singleton*.
%
%      DEERANALYSIS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DEERANALYSIS.M with the given input arguments.
%
%      DEERANALYSIS('Property','Value',...) creates a new DEERANALYSIS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before DeerAnalysis_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to DeerAnalysis_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
%
%      Modified by H. Clark Hyde, 01/07/2011, hchyde@uchicago.edu
%        *Added feature to impose distance ratio constraints in 2-gaussian or 2-rice fits
%        *New model fit constraint options are 'kmin' and 'kmax' in parameter box
%        *Fixed error conditions that occur when no data is loaded, for example:
%           -user_model_list_Callback
%           -save_Callback
%           -select_model_Callback
%        *Fixed error conditions that occur with Tikhonov Reg. par.
%         increment/decrement prior to L curve fit.
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help DeerAnalysis

% Last Modified by GUIDE v2.5 17-Jun-2018 22:09:47

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @DeerAnalysis_OpeningFcn, ...
                   'gui_OutputFcn',  @DeerAnalysis_OutputFcn, ...
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


% --- Executes just before DeerAnalysis is made visible.
function DeerAnalysis_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to DeerAnalysis (see VARARGIN)

% Choose default command line output for DeerAnalysis
handles.output = hObject;

figname='DeerAnalysis 2019 - [no data]'; % tell user, which file is current
set(handles.main_figure,'Name',figname);

whereami = which('DeerAnalysis.m');
mypath = fileparts(whereami);
tabexist = true;
if ~exist(fullfile(mypath,'pake_base40.mat'),'file')
    tabexist = false;
end
if ~exist(fullfile(mypath,'kernel64.mat'),'file')
    tabexist = false;
end
if ~exist(fullfile(mypath,'kernel128.mat'),'file')
    tabexist = false;
end
if ~exist(fullfile(mypath,'kernel256.mat'),'file')
    tabexist = false;
end
if ~exist(fullfile(mypath,'kernel512.mat'),'file')
    tabexist = false;
end
if ~exist(fullfile(mypath,'kernel1024.mat'),'file')
    tabexist = false;
end
if ~exist(fullfile(mypath,'kernel2048.mat'),'file')
    tabexist = false;
end
if ~exist(fullfile(mypath,'pake_base_tikh_512.mat'),'file')
    tabexist = false;
end
if ~tabexist
    getback = pwd;
    cd(mypath);
    make_tables;
    make_bas_Tikh;
    cd(getback);
end

% check window size and screen size and resize if required

set(handles.main_figure,'Units','pixels');
OuterPosition = get(handles.main_figure,'OuterPosition');
OuterPosition0 = OuterPosition;
fig_width = OuterPosition(3);
fig_height = OuterPosition(4);
root_units = get(0,'units');
set(0,'units','pixels');
display_size = get(0,'screensize');
set(0,'units',root_units);
disp_width = display_size(3);
disp_height = display_size(4);
size_ratio = 0.9*min([disp_width/fig_width,disp_height/fig_height]);
if size_ratio < 1
    OuterPosition(3) = floor(size_ratio*fig_width);
    OuterPosition(4) = floor(size_ratio*fig_height);
    shiftx = OuterPosition0(3)-OuterPosition(3) + floor(0.05*disp_width);
    shifty = OuterPosition0(4)-OuterPosition(4) + floor(0.05*disp_height);
    OuterPosition(1) = OuterPosition(1) - shiftx;
    OuterPosition(2) = OuterPosition(2) + shifty;
    set(handles.main_figure,'OuterPosition',OuterPosition);
end
set(handles.main_figure,'Units','characters');

handles = increaseFontSizesIfReq(handles);

handles.margin = 30; % margin for detached plots in pixels
handles.add_margin = 30; % additional margin for axes labels
handles.regpar_edit_strformat = '%0.2e'; % string format for reg.param. display

handles.time_grain = 4; % smallest allowed time increment in nanoseconds

% all data sets are empty on new start
handles.A_texp=linspace(0,1600,201); % time axis experimental data set A
handles.A_vexp=[]; % complex time-domain data set A
handles.A_vb=[];
handles.A_dipevo=[];
handles.A_cluster=[];
handles.A_bckg=[];
handles.A_tdip=[];
handles.A_r=[];
handles.A_distr=[];
handles.A_low=[];
handles.A_high=[];
handles.A_sim=[];
handles.A_spc=[];
handles.A_simspc=[];
handles.A_ny=[];
handles.A_depth=0.3;
handles.B_texp=[]; % same for data set B
handles.B_vexp=[];
handles.B_vb=[];
handles.B_dipevo=[];
handles.B_cluster=[];
handles.B_tdip=[];
handles.B_r=[];
handles.B_distr=[];
handles.B_spc=[];
handles.B_ny=[];
handles.B_sim=[];
handles.B_simspc=[];
handles.texp=[];
handles.vexp=[];
handles.t_orig=[];
handles.v_orig=[];
handles.bas_name='(none)';
handles.polynomial=[1 0]; % background polynomial
handles.dist_est=0.1;
handles.imo=0;
handles.spins_per_object=3;
handles.locked_loaded=false;
handles.Tikh_nfp_vector=[]; % Tikhonov number of effective free parameters

% Calibration values
handles.calib_density=6.46;
handles.calib_numspins=-0.5418;

% all memory flags and correction variables are reset
handles.updated=0;
handles.Tikh_updated=0;
handles.model_updated=0;
handles.imag_memory=0;
handles.phase=0;
handles.print_request=0;
handles.saved=1;
handles.validation_mode=0;
handles.theor=0;
handles.min_constraint=1.5;
handles.max_constraint=8;
handles.non_negativity = true;

imi=imread('DeerAnalysis_Logo.png','png');
axes(handles.original_data);
cla;
image(imi); 
axis off;
axis equal

% Default values internal status variables
ksize=512;
handles.kernel_size=ksize;
p_string=sprintf('%d',ksize); % display
fname=sprintf('%s%s','kernel',p_string); % generate filename
load(fname,'base','crosstalk','tnorm','ny','t');
handles.APT_kernel=base;
handles.APT_crosstalk=crosstalk;
handles.APT_norm=tnorm;
handles.APT_ny=ny;
handles.APT_t=t;
handles.zt_min_shift=1; % minimum shift in zero time (in ns)

handles.noise_averaging=false;

CData = get_detach_icon;
set(handles.pushbutton_detach_original,'CData',CData/255);
set(handles.pushbutton_detach_dipolar,'CData',CData/255);
set(handles.pushbutton_detach_distribution,'CData',CData/255);
handles.detach_icon = CData/255;
CData = get_attach_icon;
handles.attach_icon = CData/255;
handles.original_attached = true;
handles.dipolar_attached = true;
handles.distribution_attached = true;
handles.original_fig = [];
handles.dipolar_fig = [];
handles.distribution_fig = [];
handles.p1 = 1;
handles.p2 = 0;
handles.p3 = 0;
handles.bandwidth = 16;

load('pake_base40.mat','base','t','r','wd');
kernel=base-ones(size(base)); % kernel format for pcf2deer
handles.Pake_kernel=kernel;
handles.Pake_t=t;
handles.Pake_r=r;
handles.Pake_wd=wd;

load('pake_base_tikh_512','L','U','V','X','kernel','r','sm','t');
handles.Tikh_kernel = kernel';
handles.Tikh_r = r;
handles.Tikh_t = t;
handles.Tikh_U = U;
handles.Tikh_sm = sm;
handles.Tikh_X = X;
handles.Tikh_V = V;
handles.Tikh_L = L;
handles.regpar_vector = [];

rpv=load('Lcurve_abscissa.dat');
handles.rpv=rpv';
%handles.rpv=[1e-3,3.1623e-3,1e-2,3.1623e-2,1e-1,3.1623e-1,1,3.1623,1e1,3.1623e1,1e2,3.1623e2,1e3,3.1623e3,1e4,3.1623e-4,1e5]; % regularization parameters for L curve computation
handles.Lcurve_rho=zeros(1,length(rpv)); % norm of residual
handles.Lcurve_eta=zeros(1,length(rpv)); % (semi)norm of second derivative of distance distribution
handles.Lcurve_distr=zeros(length(handles.rpv),length(r));
handles.Lcurve_low=zeros(length(handles.rpv),length(r));
handles.Lcurve_high=zeros(length(handles.rpv),length(r));
handles.Lcurve_sim=[];
handles.regpar_opt_Lc=1;
handles.regpar_opt_AIC=1;
handles.regpar_opt_GCV=1;

zf=load('zero_filling.dat');
handles.zf=zf;

% Default values for Control panel: Data sets (only on program start, not
% on reload)
handles.log_Tikh=1;
set(handles.autophase,'Value',1);
set(handles.reset_on_load,'Value',1);
set(handles.format_elexsys,'Value',1);
set(handles.format_winepr,'Value',0);
set(handles.format_ascii,'Value',0);
set(handles.radiobutton_DeerAnalysis,'Value',0);
handles.ASCII_t_column=1;
pstr=sprintf('%d',handles.ASCII_t_column); % display
set(handles.ASCII_time,'String',pstr);
handles.ASCII_real_column=2;
pstr=sprintf('%d',handles.ASCII_real_column); % display
set(handles.ASCII_real,'String',pstr);
handles.ASCII_imag_column=3;
pstr=sprintf('%d',handles.ASCII_imag_column); % display
set(handles.ASCII_imaginary,'String',pstr);

handles.man_k=0.2;
handles.man_depth=0.3;

% Check for fit models

callpath=which('DeerAnalysis.m');
deer_root=callpath(1:length(callpath)-length('DeerAnalysis.m'));
models=dir(fullfile(deer_root,'models/*.m'));
handles.model_path=fullfile(deer_root,'models/');
addpath(handles.model_path);
[m,~]=size(models);
model_list=[];
if m>=1
    for k=1:m
        item=models(k).name;
        model_list=[model_list item(1:length(item)-2)];
        if k<m, model_list=[model_list '|']; end
    end
end
set(handles.user_model_list,'String',model_list);

% determine if statistics toolbox is present
% and switch off all HC Hyde features if they might be unsafe
n = 'finv'; % 'finv';
pat = '(?<=^.+[\\/]toolbox[\\/])[^\\/]+';
toolbox = regexp(which(n), pat, 'match', 'once');
vers = ver(toolbox);
if ~isempty(vers)
    set(handles.checkbox_statistics,'Enable','on');
end
fmp = which('fmincon');
if isempty(fmp)
    set(handles.checkbox_statistics,'Enable','off');
end
handles.fit_constrained = false;

% prepare for comparative mode
handles.A_curr_r = [];
handles.A_curr_distr = [];
handles.A_curr_mode = ' ';
handles.A_prev_r = [];
handles.A_prev_distr = [];
handles.A_prev_mode = ' ';
handles.new_distr = 1;

% DEERNet stuff
handles.deernet_path = fullfile(deer_root,'deernet/');
addpath(genpath(handles.deernet_path));
handles.net_set = 'net_set_any_peaks';
handles.net_set_definitions = 'https://www.dropbox.com/s/z7zzvlxnash3kn3/current_netsets.txt?dl=0';
netsets = get_netset_definitions('distributed_netsets.txt');
existing_nets = 0; 
for nd = 1:length(netsets)
    if exist(sprintf('%s%c1.mat',netsets(nd).directory,filesep),'file')
        existing_nets = existing_nets + 1;
        menu{existing_nets} = netsets(nd).name;
    end
end
% menu{existing_nets + 1} = '<check for updates>'; % needs adaptation of
% path
if existing_nets == 0
    menu{1} = '<No net sets>';
end
set(handles.popupmenu_deernet,'String',menu);
if existing_nets > 0
    set(handles.pushbutton_deernet,'String','Compute');
    set(handles.pushbutton_deernet,'TooltipString','Starts neural network analysis');
else
    set(handles.pushbutton_deernet,'String','Download');
    set(handles.pushbutton_deernet,'TooltipString','Downloads neural network sets from server');
end

handles.net_sets = netsets;

handles.new_bckg = 1;

% Remaining default behaviour of GUI
handles=set_defaults(handles);

warning off MATLAB:polyfit:RepeatedPointsOrRescale

% Include working directory into path
path(path,pwd);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes DeerAnalysis wait for user response (see UIRESUME)
% uiwait(handles.main_figure);


% --- Outputs from this function are returned to the command line.
function varargout = DeerAnalysis_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in dual_display.
function dual_display_Callback(hObject, eventdata, handles)
% hObject    handle to dual_display (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of dual_display
flag=get(hObject,'Value');
if flag
    set(handles.imaginary,'String','mod. depth scaling');
    set(handles.imaginary,'TooltipString','Automatic modulation depth scaling on/off');
    handles.imag_memory=get(handles.imaginary,'Value');
    set(handles.imaginary,'Value',0);
    handles.mask=ones(size(handles.A_r));
    set(handles.distr_suppress,'Enable','off');
else
    set(handles.imaginary,'String','imaginary');
    set(handles.imaginary,'TooltipString','Imaginary trace on/off');
    set(handles.imaginary,'Value',handles.imag_memory);
    set(handles.distr_suppress,'Enable','on');
end
% Update handles structure
guidata(hObject, handles);
update_DA(handles);

% --- Executes on button press in imaginary.
function imaginary_Callback(hObject, eventdata, handles)
% hObject    handle to imaginary (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of imaginary

update_DA(handles);

% --- Executes on button press in zt_default.
function zt_default_Callback(hObject, eventdata, handles)
% hObject    handle to zt_default (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~isfield(handles,'source_file'), set(handles.status_line,'String','### Load data file ###'); return; end
zt0=handles.zerotime;
t=handles.t_orig;
z=handles.v_orig;
phi=handles.phase*pi/180;
zt=get_zerotime(handles,t,real(z*exp(1i*phi)));
handles.zerotime=zt;
if zt~=zt0 
    handles.updated=0;     
    handles.validation_mode=0;
end
% Update handles structure
guidata(hObject, handles);
update_DA(handles);


% --- Executes on button press in zt_minus.
function zt_minus_Callback(hObject, eventdata, handles)
% hObject    handle to zt_minus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~isfield(handles,'source_file'), set(handles.status_line,'String','### Load data file ###'); return; end
[v,handles]=shift_cursor(handles,handles.zt_edit,min(handles.t_orig),max(handles.t_orig),'%d',handles.zerotime,-handles.zt_min_shift);
if v~=handles.zerotime
	handles.zerotime=v;
    handles.updated=0;
    handles.validation_mode=0;
    handles.cutoff = max(handles.t_orig)-handles.zerotime;
    set(handles.cutoff_edit,'String',sprintf('%i',handles.cutoff));
end
% Update handles structure
guidata(hObject, handles);
update_DA(handles);

% --- Executes during object creation, after setting all properties.
function zt_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to zt_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function zt_edit_Callback(hObject, eventdata, handles)
% hObject    handle to zt_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of zt_edit as text
%        str2double(get(hObject,'String')) returns contents of zt_edit as a double

if ~isfield(handles,'source_file'), set(handles.status_line,'String','### Load data file ###'); return; end
[v,handles]=edit_update(handles,hObject,min(handles.t_orig),max(handles.t_orig),0,'%d',1);
if v~=handles.zerotime
	handles.zerotime=v;
    handles.updated=0;
    handles.validation_mode=0;
    handles.cutoff = max(handles.t_orig)-handles.zerotime;
    set(handles.cutoff_edit,'String',sprintf('%i',handles.cutoff));
end
% Update handles structure
guidata(hObject, handles);
update_DA(handles);

% --- Executes on button press in zt_plus.
function zt_plus_Callback(hObject, eventdata, handles)
% hObject    handle to zt_plus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~isfield(handles,'source_file'), set(handles.status_line,'String','### Load data file ###'); return; end
[v,handles]=shift_cursor(handles,handles.zt_edit,min(handles.t_orig),max(handles.t_orig),'%d',handles.zerotime,handles.zt_min_shift);
if v~=handles.zerotime
	handles.zerotime=v;
    handles.updated=0;
    handles.validation_mode=0;
    handles.cutoff = max(handles.t_orig)-handles.zerotime;
    set(handles.cutoff_edit,'String',sprintf('%i',handles.cutoff));
end
% Update handles structure
guidata(hObject, handles);
update_DA(handles);


% --- Executes on button press in phase_default.
function phase_default_Callback(hObject, eventdata, handles)
% hObject    handle to phase_default (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~isfield(handles,'source_file'), set(handles.status_line,'String','### Load data file ###'); return; end
ph0=handles.phase;
handles=get_phase(handles);
if handles.phase~=ph0
    handles.updated=0; 
    handles.validation_mode=0;
end
% Update handles structure
guidata(hObject, handles);
update_DA(handles);

% --- Executes on button press in phase_minus.
function phase_minus_Callback(hObject, eventdata, handles)
% hObject    handle to phase_minus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~isfield(handles,'source_file'), set(handles.status_line,'String','### Load data file ###'); return; end
[v,handles]=shift_cursor(handles,handles.phase_edit,-180,180,'%5.1f',180*handles.phase/pi,-1);
if v~=handles.phase*180/pi
  	handles.phase=pi*v/180;
    handles.updated=0;
    handles.validation_mode=0;
end
% Update handles structure
guidata(hObject, handles);
update_DA(handles);

% --- Executes during object creation, after setting all properties.
function phase_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to phase_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function phase_edit_Callback(hObject, eventdata, handles)
% hObject    handle to phase_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of phase_edit as text
%        str2double(get(hObject,'String')) returns contents of phase_edit as a double
if ~isfield(handles,'source_file')
    set(hObject,'String',handles.phase)  %return to previous value
    set(handles.status_line,'String','### Load data file ###'); 
    return;
end
[v,handles]=edit_update(handles,hObject,-180,180,0,'%5.1f',0);
if v~=handles.phase*180/pi
    handles.phase=pi*v/180;
    handles.updated=0;
    handles.validation_mode=0;
end
% Update handles structure
guidata(hObject, handles);
update_DA(handles);



% --- Executes on button press in phase_plus.
function phase_plus_Callback(hObject, eventdata, handles)
% hObject    handle to phase_plus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~isfield(handles,'source_file'), set(handles.status_line,'String','### Load data file ###'); return; end
[v,handles]=shift_cursor(handles,handles.phase_edit,-180,180,'%5.1f',180*handles.phase/pi,+1);
if v~=handles.phase*180/pi
	handles.phase=pi*v/180;
    handles.updated=0;
    handles.validation_mode=0;
end
% Update handles structure
guidata(hObject, handles);
update_DA(handles);


% --- Executes on button press in bckg_default.
function bckg_default_Callback(hObject, eventdata, handles)
% hObject    handle to bckg_default (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~isfield(handles,'source_file'), set(handles.status_line,'String','### Load data file ###'); return; end
flag=get(handles.bckg_none,'Value');
if flag, return; end
bckg0=handles.bckg_start;
set(handles.main_figure,'Pointer','watch');
drawnow;
bckg=get_t_bckg(handles);
handles.bckg_start=bckg;
pstr=sprintf('%d',bckg);
set(handles.bckg_edit,'String',pstr);
set(handles.main_figure,'Pointer','arrow');
set(handles.status_line,'String','Ready.');
if handles.bckg_start~=bckg0
    handles.updated=0;
    handles.validation_mode=0;
end
% Update handles structure
handles.new_bckg = 1;
guidata(hObject, handles);
update_DA(handles);

% --- Executes on button press in bckg_minus.
function bckg_minus_Callback(hObject, eventdata, handles)
% hObject    handle to bckg_minus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~isfield(handles,'source_file'), set(handles.status_line,'String','### Load data file ###'); return; end
[v,handles]=shift_cursor(handles,handles.bckg_edit,0,handles.cutoff-handles.zerotime-handles.dt,'%d',handles.bckg_start,-handles.dt);
if v~=handles.bckg_start
    handles.bckg_start=v;
    handles.updated=0;
    handles.validation_mode=0;
end
handles.new_bckg = 1;
% Update handles structure
guidata(hObject, handles);
update_DA(handles);

% --- Executes during object creation, after setting all properties.
function bckg_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bckg_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function bckg_edit_Callback(hObject, eventdata, handles)
% hObject    handle to bckg_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of bckg_edit as text
%        str2double(get(hObject,'String')) returns contents of bckg_edit as a double
if ~isfield(handles,'source_file')
    set(hObject,'String',handles.bckg_start)  %return to previous value
    set(handles.status_line,'String','### Load data file ###');
    return
end
[v,handles]=edit_update(handles,hObject,0,handles.cutoff,0,'%d',1);
if v~=handles.bckg_start
	handles.bckg_start=v;
    handles.updated=0;
    handles.validation_mode=0;
end
handles.new_bckg = 1;
% Update handles structure
guidata(hObject, handles);
update_DA(handles);



% --- Executes on button press in bckg_plus.
function bckg_plus_Callback(hObject, eventdata, handles)
% hObject    handle to bckg_plus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~isfield(handles,'source_file'), set(handles.status_line,'String','### Load data file ###'); return; end
[v,handles]=shift_cursor(handles,handles.bckg_edit,0,handles.cutoff-handles.zerotime-handles.dt,'%d',handles.bckg_start,handles.dt);
if v~=handles.bckg_start
	handles.bckg_start=v;
    handles.updated=0;
    handles.validation_mode=0;
end
handles.new_bckg = 1;
% Update handles structure
guidata(hObject, handles);
update_DA(handles);


% --- Executes on button press in cutoff_default.
function cutoff_default_Callback(hObject, eventdata, handles)
% hObject    handle to cutoff_default (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~isfield(handles,'source_file'), set(handles.status_line,'String','### Load data file ###'); return; end
cutoff0=handles.cutoff;
handles.cutoff=handles.cutoff_suggestion;
if handles.cutoff~=cutoff0
    pstr=sprintf('%d',handles.cutoff);
    set(handles.cutoff_edit,'String',pstr);
    handles.updated=0;
    handles.validation_mode=0;
end
handles.new_bckg = 1;
% Update handles structure
guidata(hObject, handles);
update_DA(handles);

% --- Executes on button press in cutoff_minus.
function cutoff_minus_Callback(hObject, eventdata, handles)
% hObject    handle to cutoff_minus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~isfield(handles,'source_file'), set(handles.status_line,'String','### Load data file ###'); return; end
texp=pre_process(handles.t_orig,handles.v_orig,handles.phase,handles.imo,handles.zerotime,handles.dt);
[v,handles]=shift_cursor(handles,handles.cutoff_edit,handles.zerotime,max(texp),'%d',handles.cutoff,-handles.dt);
if v~=handles.cutoff
    handles.cutoff=v;
    handles.updated=0;
    handles.validation_mode=0;
end
handles.new_bckg = 1;
% Update handles structure
guidata(hObject, handles);
update_DA(handles);



% --- Executes during object creation, after setting all properties.
function cutoff_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cutoff_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function cutoff_edit_Callback(hObject, eventdata, handles)
% hObject    handle to cutoff_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of cutoff_edit as text
%        str2double(get(hObject,'String')) returns contents of cutoff_edit as a double
if ~isfield(handles,'source_file')
    set(hObject,'String',handles.cutoff)  %return to previous value
    set(handles.status_line,'String','### Load data file ###');
    return
end
texp=pre_process(handles.t_orig,handles.v_orig,handles.phase,handles.imo,handles.zerotime,handles.dt);
[v,handles]=edit_update(handles,hObject,handles.bckg_start+handles.zerotime,max(texp),max(texp),'%d',1);
if v~=handles.cutoff
	handles.cutoff=v;
    handles.updated=0;
    handles.validation_mode=0;
end
handles.new_bckg = 1;
% Update handles structure
guidata(hObject, handles);
update_DA(handles);




% --- Executes on button press in cutoff_plus.
function cutoff_plus_Callback(hObject, eventdata, handles)
% hObject    handle to cutoff_plus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~isfield(handles,'source_file'), set(handles.status_line,'String','### Load data file ###'); return; end
texp =pre_process(handles.t_orig,handles.v_orig,handles.phase,handles.imo,handles.zerotime,handles.dt);
[v,handles]=shift_cursor(handles,handles.cutoff_edit,handles.zerotime,max(texp),'%d',handles.cutoff,handles.dt);
if v~=handles.cutoff
    handles.cutoff=v;
    handles.updated=0;
    handles.validation_mode=0;
end
handles.new_bckg = 1;
% Update handles structure
guidata(hObject, handles);
update_DA(handles);

% --- Executes on button press in dip_time_domain.
function dip_time_domain_Callback(hObject, eventdata, handles)
% hObject    handle to dip_time_domain (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of dip_time_domain
flag=get(hObject,'Value');
if flag
    set(handles.dip_spectrum,'Value',0);
    set(handles.text_form_factor,'String','Form factor');
end
guidata(hObject,handles);
update_DA(handles);



% --- Executes on button press in dip_spectrum.
function dip_spectrum_Callback(hObject, eventdata, handles)
% hObject    handle to dip_spectrum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of dip_spectrum
flag=get(hObject,'Value');
if flag
    set(handles.dip_time_domain,'Value',0);
    set(handles.text_form_factor,'String','Dipolar spectrum');
end
guidata(hObject,handles);
update_DA(handles);


% --- Executes on button press in long_pass.
function long_pass_Callback(hObject, eventdata, handles)
% hObject    handle to long_pass (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of long_pass

handles.updated=0;
handles.validation_mode=0;

% Update handles structure
guidata(hObject, handles);
update_DA(handles);


% --- Executes during object creation, after setting all properties.
function long_pass_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to long_pass_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function long_pass_edit_Callback(hObject, eventdata, handles)
% hObject    handle to long_pass_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of long_pass_edit as text
%        str2double(get(hObject,'String')) returns contents of long_pass_edit as a double

v0=handles.longpass_min;
[v,handles]=edit_update(handles,hObject,1.0,2,1.5,'%5.2f',0);
handles.longpass_min=v;
if v~=v0
	handles.updated=0;
    handles.validation_mode=0;
end
% Update handles structure
guidata(hObject, handles);
update_DA(handles);


% --- Executes on button press in exci_bandwidth_corr.
function exci_bandwidth_corr_Callback(hObject, eventdata, handles)
% hObject    handle to exci_bandwidth_corr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of exci_bandwidth_corr

handles.updated=0;
handles.validation_mode=0;
% Update handles structure
guidata(hObject, handles);
update_DA(handles);

% --- Executes during object creation, after setting all properties.
function exci_bandwidth_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to exci_bandwidth_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function exci_bandwidth_edit_Callback(hObject, eventdata, handles)
% hObject    handle to exci_bandwidth_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of exci_bandwidth_edit as text
%        str2double(get(hObject,'String')) returns contents of exci_bandwidth_edit as a double

v0=handles.bandwidth;
[v,handles]=edit_update(handles,hObject,1.0,1e6,16,'%0.3g',0);
handles.bandwidth=v;
if v~=v0
    handles.updated=0;
    handles.validation_mode=0;
end
% Update handles structure
guidata(hObject, handles);
update_DA(handles);

% --- Executes on button press in dip_expand.
function dip_expand_Callback(hObject, eventdata, handles)
% hObject    handle to dip_expand (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.zoom=handles.zoom*sqrt(2);
if handles.zoom>1000 
    handles.zoom=1000;
end
pstr=sprintf('%5.1f',handles.zoom);
set(handles.zoom_edit,'String',pstr);
% Update handles structure
guidata(hObject, handles);
update_DA(handles);



% --- Executes on button press in dip_full.
function dip_full_Callback(hObject, eventdata, handles)
% hObject    handle to dip_full (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.zoom=1;
set(handles.zoom_edit,'String','1.0');
% Update handles structure
guidata(hObject, handles);
update_DA(handles);


% --- Executes on button press in dip_shrink.
function dip_shrink_Callback(hObject, eventdata, handles)
% hObject    handle to dip_shrink (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.zoom=handles.zoom/sqrt(2);
if handles.zoom<1
    handles.zoom=1;
end
pstr=sprintf('%5.1f',handles.zoom);
set(handles.zoom_edit,'String',pstr);
% Update handles structure
guidata(hObject, handles);
update_DA(handles);


% --- Executes during object creation, after setting all properties.
function zoom_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to zoom_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function zoom_edit_Callback(hObject, eventdata, handles)
% hObject    handle to zoom_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of zoom_edit as text
%        str2double(get(hObject,'String')) returns contents of zoom_edit as a double
[v,handles]=edit_update(handles,hObject,1.0,1000,1.0,'%5.1f',0);
handles.zoom=v;
% Update handles structure
guidata(hObject, handles);
update_DA(handles);


% --- Executes on button press in load.
function load_Callback(hObject, eventdata, handles)
% hObject    handle to load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handback = get_dataset(handles);
set(handles.distr_suppress,'Enable','on');
guidata(hObject,handback);
update_DA(handback);


% --- Executes on button press in load_series.
function load_series_Callback(hObject, eventdata, handles)
% hObject    handle to load_series (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles=process_list(handles);
guidata(hObject,handles);
update_DA(handles);


% --- Executes on button press in format_elexsys.
function format_elexsys_Callback(hObject, eventdata, handles)
% hObject    handle to format_elexsys (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of format_elexsys
flag=get(hObject,'Value');
if flag
    set(handles.format_winepr,'Value',0);
    set(handles.format_ascii,'Value',0);
    set(handles.radiobutton_DeerAnalysis,'Value',0);
    set(handles.checkbox_no_analysis,'Enable','off');
    set(handles.checkbox_no_analysis,'Value',0);
end
guidata(hObject,handles);


% --- Executes on button press in format_winepr.
function format_winepr_Callback(hObject, eventdata, handles)
% hObject    handle to format_winepr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of format_winepr
flag=get(hObject,'Value');
if flag
    set(handles.format_elexsys,'Value',0);
    set(handles.format_ascii,'Value',0);
    set(handles.radiobutton_DeerAnalysis,'Value',0);
    set(handles.checkbox_no_analysis,'Enable','off');
    set(handles.checkbox_no_analysis,'Value',0);
end
guidata(hObject,handles);

% --- Executes on button press in format_ascii.
function format_ascii_Callback(hObject, eventdata, handles)
% hObject    handle to format_ascii (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of format_ascii
flag=get(hObject,'Value');
if flag
    set(handles.format_elexsys,'Value',0);
    set(handles.format_winepr,'Value',0);
    set(handles.radiobutton_DeerAnalysis,'Value',0);
    set(handles.checkbox_no_analysis,'Enable','off');
    set(handles.checkbox_no_analysis,'Value',0);
end
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function edit8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function edit8_Callback(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit8 as text
%        str2double(get(hObject,'String')) returns contents of edit8 as a double


% --- Executes during object creation, after setting all properties.
function ASCII_time_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ASCII_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function ASCII_time_Callback(hObject, eventdata, handles)
% hObject    handle to ASCII_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ASCII_time as text
%        str2double(get(hObject,'String')) returns contents of ASCII_time as a double

[v,handles]=edit_update(handles,hObject,1,100,1,'%d',1);
handles.ASCII_t_column=v;
% Update handles structure
guidata(hObject, handles);
update_DA(handles);

% --- Executes during object creation, after setting all properties.
function ASCII_real_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ASCII_real (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function ASCII_real_Callback(hObject, eventdata, handles)
% hObject    handle to ASCII_real (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ASCII_real as text
%        str2double(get(hObject,'String')) returns contents of ASCII_real as a double

[v,handles]=edit_update(handles,hObject,1,100,2,'%d',1);
handles.ASCII_real_column=v;
% Update handles structure
guidata(hObject, handles);
update_DA(handles);



% --- Executes during object creation, after setting all properties.
function ASCII_imaginary_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ASCII_imaginary (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function ASCII_imaginary_Callback(hObject, eventdata, handles)
% hObject    handle to ASCII_imaginary (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ASCII_imaginary as text
%        str2double(get(hObject,'String')) returns contents of ASCII_imaginary as a double

[v,handles]=edit_update(handles,hObject,1,100,1,'%d',1);
handles.ASCII_imag_column=v;
% Update handles structure
guidata(hObject, handles);
update_DA(handles);



% --- Executes on button press in save.
function save_Callback(hObject, eventdata, handles)
% hObject    handle to save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

return_path=pwd;
cd(handles.project_dir);
suggestion=[handles.bas_name '_res.txt']; 
[fname,pname]=uiputfile('*.txt','Save results',suggestion);
% Remove (last) extension, if any
s=strfind(fname,'.');
if ~isempty(s)
	fname=fname(1:s(length(s))-1);
end
% Remove suffix '_res', if present
s=strfind(fname,'_res');
if ~isempty(s)
	fname=fname(1:s(length(s))-1);
end
handles=save_result(handles,fname,pname);
cd(return_path);
guidata(hObject, handles);

% --- Executes on button press in Copy.
function Copy_Callback(hObject, eventdata, handles)
% hObject    handle to Copy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.print_request=1;
% Update handles structure
guidata(hObject, handles);
update_DA(handles);


% --- Executes on button press in Close.
function Close_Callback(hObject, eventdata, handles)
% hObject    handle to Close (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

main_figure_CloseRequestFcn(handles.main_figure,eventdata,handles);

% shh = get(0,'ShowHiddenHandles');
% set(0,'ShowHiddenHandles','on');
% delete(findobj(get(0,'Children'),'flat','Tag','StatusBar'));
% set(0,'ShowHiddenHandles',shh);
% 
% if handles.saved,
% 	close(handles.main_figure);
% else
%    answer=questdlg('Do you want to save now?','Fit data not saved');
%    if strncmp(answer,'Yes',2);
% 		return_path=pwd;
% 		cd(handles.project_dir);
% 		suggestion=[handles.bas_name '_res.txt']; 
% 		[fname,pname]=uiputfile('*.txt','Save results',suggestion);
% 		% Remove (last) extension, if any
% 		s=strfind(fname,'.');
% 		if ~isempty(s),
% 			fname=fname(1:s(length(s))-1);
% 		end;
% 		% Remove suffix '_res', if present
% 		s=strfind(fname,'_res');
% 		if ~isempty(s),
% 			fname=fname(1:s(length(s))-1);
% 		end;
% 		handles=save_result(handles,fname,pname);
% 		cd(return_path);
%     	close(handles.main_figure);
%    end;
%    if strncmp(answer,'No',2);
%        close(handles.main_figure);
%    end;
% end;

% --- Executes on button press in L_curve.
function L_curve_Callback(hObject, eventdata, handles)
% hObject    handle to L_curve (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of L_curve
% Update handles structure
guidata(hObject, handles);
update_DA(handles);


% --- Executes on button press in regpar_default_Lcorner.
function regpar_default_Lcorner_Callback(hObject, eventdata, handles)
% hObject    handle to regpar_default_Lcorner (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~isfield(handles,'source_file'), set(handles.status_line,'String','### Load data file ###'); return; end
if isempty(handles.Lcurve_sim), set(handles.status_line,'String','### Tikhonov fit required. ###'); return; end
if numel(handles.A_r)~=numel(handles.A_distr), set(handles.status_line,'String','### Tikhonov L curve fit required. ###'); return; end
if numel(handles.regpar_vector)==1, set(handles.status_line,'String','### Tikhonov L curve fit required. ###'); return; end

idxLc = handles.regpar_opt_Lc;
handles.regpar_sel = idxLc;
handles.regpar = handles.regpar_vector(idxLc);

set(handles.regpar_edit,'String',num2str(handles.regpar,handles.regpar_edit_strformat));
set(handles.status_line,'String','Recomputing...');
set(handles.main_figure,'Pointer','watch');
drawnow

[r,distr] = get_Tikhonov_new(handles,handles.regpar);

set(handles.status_line,'String','Simulating form factor...');
% check if excitation bandwidth correction is selected
exBWcorr = get(handles.exci_bandwidth_corr,'Value'); 
if exBWcorr
    [sim,scale] = deer_sim(r,distr,handles.A_tdip,handles.bandwidth);
elseif length(handles.A_tdip) > 1024
    sim = deer_sim(r,distr,handles.A_tdip);
    scale = 1;
else
    sim = get_td_fit(handles,r,distr);
    scale = 1;
end
handles.moddepth_suppression = scale;
handles.A_sim = sim;

set(handles.status_line,'String','Ready.');
set(handles.main_figure,'Pointer','arrow');
drawnow
handles.A_r=r;
handles.A_distr=distr';
handles.A_low=distr';
handles.A_high=distr';
handles.mask=ones(size(handles.A_distr));
% Update handles structure
guidata(hObject, handles);
update_DA(handles);

% --- Executes on button press in regpar_AIC.
function regpar_AIC_Callback(hObject, eventdata, handles)
% hObject    handle to regpar_AIC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~isfield(handles,'source_file'), set(handles.status_line,'String','### Load data file ###'); return; end
if isempty(handles.Lcurve_sim), set(handles.status_line,'String','### Tikhonov fit required. ###'); return; end
if numel(handles.A_r)~=numel(handles.A_distr), set(handles.status_line,'String','### Tikhonov L curve fit required. ###'); return; end
if numel(handles.regpar_vector)==1, set(handles.status_line,'String','### Tikhonov L curve fit required. ###'); return; end

idxAIC = handles.regpar_opt_AIC;
handles.regpar_sel = idxAIC;
handles.regpar = handles.regpar_vector(idxAIC);

set(handles.regpar_edit,'String',num2str(handles.regpar,handles.regpar_edit_strformat));
set(handles.status_line,'String','Recomputing...');
set(handles.main_figure,'Pointer','watch');
drawnow

[r,distr] = get_Tikhonov_new(handles,handles.regpar);

set(handles.status_line,'String','Simulating form factor...');
% check if excitation bandwidth correction is selected
exBWcorr = get(handles.exci_bandwidth_corr,'Value'); 
if exBWcorr
    [sim,scale] = deer_sim(r,distr,handles.A_tdip,handles.bandwidth);
elseif length(handles.A_tdip) > 1024
    sim = deer_sim(r,distr,handles.A_tdip);
    scale = 1;
else
    sim = get_td_fit(handles,r,distr);
    scale = 1;
end
handles.moddepth_suppression = scale;
handles.A_sim = sim;

set(handles.status_line,'String','Ready.');
set(handles.main_figure,'Pointer','arrow');
drawnow
handles.A_r=r;
handles.A_distr=distr';
handles.A_low=distr';
handles.A_high=distr';
handles.mask=ones(size(handles.A_distr));
% Update handles structure
guidata(hObject, handles);
update_DA(handles);

% --- Executes on button press in regpar_up.
function regpar_up_Callback(hObject, eventdata, handles)
% hObject    handle to regpar_up (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~isfield(handles,'source_file'), set(handles.status_line,'String','### Load data file ###'); return; end
if numel(handles.A_r)~=numel(handles.A_distr), set(handles.status_line,'String','### Tikhonov L curve fit required. ###'); return; end
if length(handles.regpar_vector)<=1, set(handles.status_line,'String','### Tikhonov L curve fit required. ###'); return; end
v = max(handles.regpar_sel-1,1);
handles.regpar_sel=v;
handles.regpar=handles.regpar_vector(v);

set(handles.regpar_edit,'String',num2str(handles.regpar,handles.regpar_edit_strformat));
set(handles.status_line,'String','Recomputing...');
set(handles.main_figure,'Pointer','watch');
drawnow
[r,distr] = get_Tikhonov_new(handles,handles.regpar);
set(handles.status_line,'String','Simulating form factor...');
exflag=get(handles.exci_bandwidth_corr,'Value'); % check, if excitation bandwidth correction is selected
if exflag
    [sim,sc]=deer_sim(r,distr,handles.A_tdip,handles.bandwidth);
elseif length(handles.A_tdip) > 1024
    sim = deer_sim(r,distr,handles.A_tdip);
    sc = 1;
else
    sim=get_td_fit(handles,r,distr);
    sc = 1;
end
handles.moddepth_suppression=sc;
handles.A_sim=sim;
set(handles.status_line,'String','Ready.');
set(handles.main_figure,'Pointer','arrow');
drawnow
handles.A_r=r;
handles.A_distr=distr';
handles.A_low=distr';
handles.A_high=distr';
handles.mask=ones(size(handles.A_distr));
% Update handles structure
guidata(hObject, handles);
update_DA(handles);



% --- Executes during object creation, after setting all properties.
function regpar_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to regpar_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function regpar_edit_Callback(hObject, eventdata, handles)
% hObject    handle to regpar_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of regpar_edit as text
%        str2double(get(hObject,'String')) returns contents of regpar_edit as a double

if ~isfield(handles,'source_file')
    set(hObject,'String',handles.regpar)  %return to previous value
    set(handles.status_line,'String','### Load data file ###');
    return;
end
v0=handles.regpar;
[v,handles]=edit_update(handles,hObject,1e-9,1e9,1,'%0.5g',0);
handles.regpar=v;
for k=1:length(handles.regpar_vector)
    if abs(handles.regpar_vector(k)-v)/handles.regpar_vector(k)<1.0e-3
        handles.regpar=handles.regpar_vector(k);
        handles.regpar_sel = k;
        v0=v;
    end
end
if v~=v0
    set(handles.L_curve,'Value',0);
    set(handles.L_curve,'Enable','off');
    set(handles.select_L_curve,'Value',0);
end
set(handles.regpar_edit,'String',num2str(handles.regpar,handles.regpar_edit_strformat));
set(handles.status_line,'String','Recomputing...');
set(handles.main_figure,'Pointer','watch');
drawnow
[r,distr] = get_Tikhonov_new(handles,handles.regpar);
set(handles.status_line,'String','Simulating form factor...');
exflag=get(handles.exci_bandwidth_corr,'Value'); % check, if excitation bandwidth correction is selected
if exflag
    [sim,sc]=deer_sim(r,distr,handles.A_tdip,handles.bandwidth);
elseif length(handles.A_tdip) > 1024
    sim = deer_sim(r,distr,handles.A_tdip);
    sc = 1;
else
    sim=get_td_fit(handles,r,distr);
    sc = 1;
end
handles.moddepth_suppression=sc;
handles.A_sim=sim;
set(handles.status_line,'String','Ready.');
set(handles.main_figure,'Pointer','arrow');
drawnow
handles.A_r=r;
handles.A_distr=distr';
handles.A_low=distr';
handles.A_high=distr';
handles.mask=ones(size(handles.A_distr));
% Update handles structure
guidata(hObject, handles);
update_DA(handles);

% --- Executes on button press in regpar_down.
function regpar_down_Callback(hObject, eventdata, handles)
% hObject    handle to regpar_down (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~isfield(handles,'source_file'), set(handles.status_line,'String','### Load data file ###'); return; end
if numel(handles.A_r)~=numel(handles.A_distr), set(handles.status_line,'String','### Tikhonov L curve fit required. ###'); return; end
if length(handles.regpar_vector)==1, set(handles.status_line,'String','### Tikhonov L curve fit required. ###'); return; end
m=length(handles.regpar_vector);
v=min(handles.regpar_sel+1,m);
handles.regpar_sel=v;
handles.regpar=handles.regpar_vector(v);
set(handles.regpar_edit,'String',num2str(handles.regpar,handles.regpar_edit_strformat));
set(handles.status_line,'String','Recomputing...');
set(handles.main_figure,'Pointer','watch');
drawnow
[r,distr] = get_Tikhonov_new(handles,handles.regpar);
set(handles.status_line,'String','Simulating form factor...');
exflag=get(handles.exci_bandwidth_corr,'Value'); % check, if excitation bandwidth correction is selected
if exflag
    [sim,sc]=deer_sim(r,distr,handles.A_tdip,handles.bandwidth);
elseif length(handles.A_tdip) > 1024
    sim = deer_sim(r,distr,handles.A_tdip);
    sc = 1;
else
    sim=get_td_fit(handles,r,distr);
    sc = 1;
end
handles.moddepth_suppression=sc;
handles.A_sim=sim;
set(handles.status_line,'String','Ready.');
set(handles.main_figure,'Pointer','arrow');
drawnow
handles.A_r=r;
handles.A_distr=distr';
handles.A_low=distr';
handles.A_high=distr';
handles.mask=ones(size(handles.A_distr));
% Update handles structure
guidata(hObject, handles);
update_DA(handles);

% --- Executes on button press in distr_start_default.
function distr_start_default_Callback(hObject, eventdata, handles)
% hObject    handle to distr_start_default (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.rmin=1.5;
set(handles.distr_start_edit,'String','1.50');
% Update handles structure
guidata(hObject, handles);
update_DA(handles);

% --- Executes on button press in distr_start_default_minus.
function distr_start_minus_Callback(hObject, eventdata, handles)
% hObject    handle to distr_start_default_minus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[v,handles]=shift_cursor(handles,handles.distr_start_edit,0.5,handles.rmax-0.05,'%5.2f',handles.rmin,-0.05);
handles.rmin=v;
% Update handles structure
guidata(hObject, handles);
update_DA(handles);


% --- Executes during object creation, after setting all properties.
function distr_start_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to distr_start_default_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function distr_start_edit_Callback(hObject, eventdata, handles)
% hObject    handle to distr_start_default_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of distr_start_default_edit as text
%        str2double(get(hObject,'String')) returns contents of distr_start_default_edit as a double
[v,handles]=edit_update(handles,hObject,0.5,handles.rmax-0.05,1.5,'%5.2f',0);
handles.rmin=v;
% Update handles structure
guidata(hObject, handles);
update_DA(handles);

% --- Executes on button press in distr_start_default_plus.
function distr_start_plus_Callback(hObject, eventdata, handles)
% hObject    handle to distr_start_default_plus (see GCBO)
% eventdata  reserved - to be dened in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[v,handles]=shift_cursor(handles,handles.distr_start_edit,0.5,handles.rmax-0.05,'%5.2f',handles.rmin,0.05);
handles.rmin=v;
% Update handles structure
guidata(hObject, handles);
update_DA(handles);


% --- Executes on button press in distr_end_plus.
function distr_end_plus_Callback(hObject, eventdata, handles)
% hObject    handle to distr_end_plus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[v,handles]=shift_cursor(handles,handles.distr_end_edit,handles.rmin+0.05,20,'%5.2f',handles.rmax,+0.05);
handles.rmax=v;
% Update handles structure
guidata(hObject, handles);
update_DA(handles);

% --- Executes during object creation, after setting all properties.
function distr_end_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to distr_end_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function distr_end_edit_Callback(hObject, eventdata, handles)
% hObject    handle to distr_end_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of distr_end_edit as text
%        str2double(get(hObject,'String')) returns contents of distr_end_edit as a double
[v,handles]=edit_update(handles,hObject,handles.rmin+0.05,20,1.5,'%5.2f',0);
handles.rmax=v;
% Update handles structure
guidata(hObject, handles);
update_DA(handles);



% --- Executes on button press in distr_end_minus.
function distr_end_minus_Callback(hObject, eventdata, handles)
% hObject    handle to distr_end_minus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[v,handles]=shift_cursor(handles,handles.distr_end_edit,handles.rmin+0.05,20,'%5.2f',handles.rmax,-0.05);
handles.rmax=v;
% Update handles structure
guidata(hObject, handles);
update_DA(handles);


% --- Executes on button press in distr_end_default.
function distr_end_default_Callback(hObject, eventdata, handles)
% hObject    handle to distr_end_default (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.rmax=8.0;
set(handles.distr_end_edit,'String','8.00');
% Update handles structure
guidata(hObject, handles);
update_DA(handles);


% --- Executes on button press in distr_expand.
function distr_expand_Callback(hObject, eventdata, handles)
% hObject    handle to distr_expand (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Update handles structure
guidata(hObject, handles);
update_DA(handles);

% --- Executes on button press in distr_suppress.
function distr_suppress_Callback(hObject, eventdata, handles)
% hObject    handle to distr_suppress (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~isfield(handles,'source_file'), set(handles.status_line,'String','### Load data file ###'); return; end
rmi=handles.A_r-handles.rmin*ones(size(handles.A_r));
[mm,pmi]=min(abs(rmi));
rma=handles.A_r-handles.rmax*ones(size(handles.A_r));
[mm,pma]=min(abs(rma));
for k=pmi:pma
    handles.mask(k)=0;
end
mask_distr=handles.A_distr.*handles.mask;
exflag=get(handles.exci_bandwidth_corr,'Value');
APTflag=get(handles.select_APT,'Value');
set(handles.status_line,'String','Simulating DEER data...');
set(handles.main_figure,'Pointer','watch');
if exflag && ~APTflag
    sim = deer_sim(handles.A_r,mask_distr,handles.A_tdip,handles.bandwidth);
else
    sim=get_td_fit(handles,handles.A_r,mask_distr);
end
set(handles.status_line,'String','Ready.');
set(handles.main_figure,'Pointer','arrow');
modsim=ones(size(sim))-sim;
modexp=ones(size(handles.A_cluster))-handles.A_cluster;
sc=sum(modexp.*modexp)/sum(modsim.*modexp);
sim=ones(size(modsim))-sc*modsim;
handles.mask_sim=sim;
% Update handles structure
guidata(hObject, handles);
update_DA(handles);



% --- Executes on button press in autophase.
function autophase_Callback(hObject, eventdata, handles)
% hObject    handle to autophase (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of autophase


% --- Executes on button press in numspin_calibrate.
function numspin_calibrate_Callback(hObject, eventdata, handles)
% hObject    handle to numspin_calibrate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function bckg_density_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bckg_density (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function bckg_density_Callback(hObject, eventdata, handles)
% hObject    handle to bckg_density (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of bckg_density as text
%        str2double(get(hObject,'String')) returns contents of bckg_density as a double

[v,handles]=edit_update(handles,hObject,1e-6,1e6,handles.bckg_dens*handles.calib_density,'%5.2f',0);
handles.calib_density=v/handles.bckg_dens;
set(hObject,'ForegroundColor','g');
handles.new_bckg = 1;
guidata(hObject,handles);
update_DA(handles);


% --- Executes on button press in reset_on_load.
function reset_on_load_Callback(hObject, eventdata, handles)
% hObject    handle to reset_on_load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of reset_on_load


% --- Executes on button press in bckg_calibrate.
function bckg_calibrate_Callback(hObject, eventdata, handles)
% hObject    handle to bckg_calibrate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.bckg_density,'ForegroundColor','g');
guidata(hObject,handles);
update_DA(handles);

% --- Executes on button press in select_APT.
function select_APT_Callback(hObject, eventdata, handles)
% hObject    handle to select_APT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of select_APT
flag=get(hObject,'Value');
if flag
    set(handles.select_Tikhonov,'Value',0);
    set(handles.select_model,'Value',0);
    set(handles.select_deernet,'Value',0);
    handles.model_updated=1;
    set(handles.L_curve,'Value',0);
    set(handles.L_curve,'Enable','off');
    handles.new_distr = 1;
    handles.validation_mode=0;
    set(handles.error_estimate,'Value',0);
end
guidata(hObject,handles);
update_DA(handles);


% --- Executes during object creation, after setting all properties.
function DDS_filter_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DDS_filter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function DDS_filter_Callback(hObject, eventdata, handles)
% hObject    handle to DDS_filter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of DDS_filter as text
%        str2double(get(hObject,'String')) returns contents of DDS_filter as a double
[v,handles]=edit_update(handles,hObject,0.05,10,0.2,'%5.2f',0);
handles.DDS=v;
% Update handles structure
guidata(hObject, handles);
update_DA(handles);

% --- Executes on button press in select_Tikhonov.
function select_Tikhonov_Callback(hObject, eventdata, handles)
% hObject    handle to select_Tikhonov (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of select_Tikhonov
flag=get(hObject,'Value');
if flag
    handles.updated=0;
    handles.validation_mode=0;
    set(handles.select_APT,'Value',0);
    set(handles.select_model,'Value',0);
    set(handles.select_deernet,'Value',0);
    handles.model_updated=1;
    handles.new_distr = 1;
    set(handles.error_estimate,'Value',0);
end
guidata(hObject,handles);
update_DA(handles);


% --- Executes on button press in select_L_curve.
function select_L_curve_Callback(hObject, eventdata, handles)
% hObject    handle to select_L_curve (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of select_L_curve


% --- Executes on button press in select_model.
function select_model_Callback(hObject, eventdata, handles)
% hObject    handle to select_model (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of select_model
flag=get(hObject,'Value');
if flag
    handles.updated=0;
    handles.validation_mode=0;
    set(handles.select_APT,'Value',0);
    set(handles.select_Tikhonov,'Value',0);
    set(handles.select_deernet,'Value',0);
    set(handles.L_curve,'Value',0);
    set(handles.L_curve,'Enable','off');
    if isfield(handles,'source_file')
        handles=sim_user_model(handles);
    end
    handles.new_distr = 1;
    set(handles.error_estimate,'Value',0);
else
    handles.model_updated=1;
end
guidata(hObject,handles);
update_DA(handles);


% --- Executes during object creation, after setting all properties.
function user_model_list_CreateFcn(hObject, eventdata, handles)
% hObject    handle to user_model_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on selection change in user_model_list.
function user_model_list_Callback(hObject, eventdata, handles)
% hObject    handle to user_model_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns user_model_list contents as cell array
%        contents{get(hObject,'Value')} returns selected item from user_model_list

handles=read_user_model(handles);
handles.updated=0;
handles.validation_mode=0;
if isfield(handles,'source_file')
    handles.model_updated=0;
    handles=sim_user_model(handles);
    guidata(hObject,handles);
    update_DA(handles);
end


% --- Executes during object creation, after setting all properties.
function par1_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to par1_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function par1_edit_Callback(hObject, eventdata, handles)
% hObject    handle to par1_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of par1_edit as text
%        str2double(get(hObject,'String')) returns contents of par1_edit as a double

[v,handles]=edit_update(handles,hObject,handles.model_lower_bounds(1),handles.model_upper_bounds(1),handles.model_defaults(1),'%7.3f',0);
handles.model_pars(1)=v;
handles.updated=0;
handles.validation_mode=0;
if isfield(handles,'source_file')
   handles=sim_user_model(handles); 
   handles.model_updated=0;
end
% Update handles structure
guidata(hObject, handles);
update_DA(handles);


% --- Executes on button press in sel_par2.
function sel_par2_Callback(hObject, eventdata, handles)
% hObject    handle to sel_par2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of sel_par2


% --- Executes during object creation, after setting all properties.
function par2_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to par2_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function par2_edit_Callback(hObject, eventdata, handles)
% hObject    handle to par2_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of par2_edit as text
%        str2double(get(hObject,'String')) returns contents of par2_edit as a double

[v,handles]=edit_update(handles,hObject,handles.model_lower_bounds(2),handles.model_upper_bounds(2),handles.model_defaults(2),'%7.3f',0);
handles.model_pars(2)=v;
handles.updated=0;
handles.validation_mode=0;
if isfield(handles,'source_file')
   handles=sim_user_model(handles); 
end
handles.model_updated=0;
% Update handles structure
guidata(hObject, handles);
update_DA(handles);



% --- Executes on button press in sel_par3.
function sel_par3_Callback(hObject, eventdata, handles)
% hObject    handle to sel_par3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of sel_par3


% --- Executes during object creation, after setting all properties.
function par3_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to par3_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function par3_edit_Callback(hObject, eventdata, handles)
% hObject    handle to par3_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of par3_edit as text
%        str2double(get(hObject,'String')) returns contents of par3_edit as a double

[v,handles]=edit_update(handles,hObject,handles.model_lower_bounds(3),handles.model_upper_bounds(3),handles.model_defaults(3),'%7.3f',0);
handles.model_pars(3)=v;
handles.updated=0;
handles.validation_mode=0;
if isfield(handles,'source_file')
   handles=sim_user_model(handles); 
   handles.model_updated=0;
end
% Update handles structure
guidata(hObject, handles);
update_DA(handles);


% --- Executes on button press in sel_par4.
function sel_par4_Callback(hObject, eventdata, handles)
% hObject    handle to sel_par4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of sel_par4


% --- Executes during object creation, after setting all properties.
function par4_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to par4_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function par4_edit_Callback(hObject, eventdata, handles)
% hObject    handle to par4_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of par4_edit as text
%        str2double(get(hObject,'String')) returns contents of par4_edit as a double

[v,handles]=edit_update(handles,hObject,handles.model_lower_bounds(4),handles.model_upper_bounds(4),handles.model_defaults(4),'%7.3f',0);
handles.model_pars(4)=v;
handles.updated=0;
handles.validation_mode=0;
if isfield(handles,'source_file')
   handles=sim_user_model(handles); 
   handles.model_updated=0;
end
% Update handles structure
guidata(hObject, handles);
update_DA(handles);


% --- Executes on button press in sel_par5.
function sel_par5_Callback(hObject, eventdata, handles)
% hObject    handle to sel_par5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of sel_par5


% --- Executes during object creation, after setting all properties.
function par5_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to par5_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function par5_edit_Callback(hObject, eventdata, handles)
% hObject    handle to par5_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of par5_edit as text
%        str2double(get(hObject,'String')) returns contents of par5_edit as a double

[v,handles]=edit_update(handles,hObject,handles.model_lower_bounds(5),handles.model_upper_bounds(5),handles.model_defaults(5),'%7.3f',0);
handles.model_pars(5)=v;
handles.updated=0;
handles.validation_mode=0;
if isfield(handles,'source_file')
   handles=sim_user_model(handles); 
    handles.model_updated=0;
end
% Update handles structure
guidata(hObject, handles);
update_DA(handles);


% --- Executes on button press in sel_par6.
function sel_par6_Callback(hObject, eventdata, handles)
% hObject    handle to sel_par6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of sel_par6


% --- Executes during object creation, after setting all properties.
function par6_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to par6_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function par6_edit_Callback(hObject, eventdata, handles)
% hObject    handle to par6_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of par6_edit as text
%        str2double(get(hObject,'String')) returns contents of par6_edit as a double

[v,handles]=edit_update(handles,hObject,handles.model_lower_bounds(6),handles.model_upper_bounds(6),handles.model_defaults(6),'%7.3f',0);
handles.model_pars(6)=v;
handles.updated=0;
handles.validation_mode=0;
if isfield(handles,'source_file')
   handles=sim_user_model(handles); 
    handles.model_updated=0;
end
% Update handles structure
guidata(hObject, handles);
update_DA(handles);

% --- Executes on button press in sel_par7.
function sel_par7_Callback(hObject, eventdata, handles)
% hObject    handle to sel_par7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of sel_par7



function par7_edit_Callback(hObject, eventdata, handles)
% hObject    handle to par7_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of par7_edit as text
%        str2double(get(hObject,'String')) returns contents of par7_edit as a double
[v,handles]=edit_update(handles,hObject,handles.model_lower_bounds(7),handles.model_upper_bounds(7),handles.model_defaults(7),'%7.3f',0);
handles.model_pars(7)=v;
handles.updated=0;
handles.validation_mode=0;
if isfield(handles,'source_file')
   handles=sim_user_model(handles); 
    handles.model_updated=0;
end
% Update handles structure
guidata(hObject, handles);
update_DA(handles);


% --- Executes during object creation, after setting all properties.
function par7_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to par7_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in sel_par8.
function sel_par8_Callback(hObject, eventdata, handles)
% hObject    handle to sel_par8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of sel_par8



function par8_edit_Callback(hObject, eventdata, handles)
% hObject    handle to par8_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of par8_edit as text
%        str2double(get(hObject,'String')) returns contents of par8_edit as a double
[v,handles]=edit_update(handles,hObject,handles.model_lower_bounds(8),handles.model_upper_bounds(8),handles.model_defaults(8),'%7.3f',0);
handles.model_pars(8)=v;
handles.updated=0;
handles.validation_mode=0;
if isfield(handles,'source_file')
   handles=sim_user_model(handles); 
    handles.model_updated=0;
end
% Update handles structure
guidata(hObject, handles);
update_DA(handles);


% --- Executes during object creation, after setting all properties.
function par8_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to par8_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in bckg_none.
function bckg_none_Callback(hObject, eventdata, handles)
% hObject    handle to bckg_none (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of bckg_none
flag=get(hObject,'Value');
if flag
    set(handles.bckg_homogeneous,'Value',0);
    set(handles.bckg_poly,'Value',0);
    set(handles.bckg_exp,'Value',0);
    set(handles.radiobutton_deernet_bckg,'Value',0);
end
handles.updated=0;
handles.model_updated=0;
handles.validation_mode=0;
handles.new_bckg = 1;
guidata(hObject,handles);
update_DA(handles);


% --- Executes on button press in bckg_homogeneous.
function bckg_homogeneous_Callback(hObject, eventdata, handles)
% hObject    handle to bckg_homogeneous (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of bckg_homogeneous
flag=get(hObject,'Value');
if flag
    set(handles.bckg_none,'Value',0);
    set(handles.bckg_poly,'Value',0);
    set(handles.bckg_exp,'Value',0);
    set(handles.radiobutton_deernet_bckg,'Value',0);
end
handles.updated=0;
handles.model_updated=0;
handles.validation_mode=0;
handles.new_bckg = 1;
guidata(hObject,handles);
update_DA(handles);


% --- Executes during object creation, after setting all properties.
function bckg_dim_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bckg_dim_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function bckg_dim_edit_Callback(hObject, eventdata, handles)
% hObject    handle to bckg_dim_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of bckg_dim_edit as text
%        str2double(get(hObject,'String')) returns contents of bckg_dim_edit as a double
[v,handles]=edit_update(handles,hObject,0.01,20,3,'%5.2f',0);
if v~=handles.hom_dim
	handles.hom_dim=v;
    handles.updated=0;
    handles.validation_mode=0;
end
handles.new_bckg = 1;
% Update handles structure
guidata(hObject, handles);
update_DA(handles);

% --- Executes on button press in bckg_fit_dim.
function bckg_fit_dim_Callback(hObject, eventdata, handles)
% hObject    handle to bckg_fit_dim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of bckg_fit_dim

handles.updated=0;
handles.validation_mode=0;
handles.new_bckg = 1;
% Update handles structure
guidata(hObject, handles);
update_DA(handles);

% --- Executes on button press in bckg_poly.
function bckg_poly_Callback(hObject, eventdata, handles)
% hObject    handle to bckg_poly (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of bckg_poly
flag=get(hObject,'Value');
if flag
    set(handles.bckg_homogeneous,'Value',0);
    set(handles.bckg_none,'Value',0);
    set(handles.bckg_exp,'Value',0);
    set(handles.radiobutton_deernet_bckg,'Value',0);
end
handles.updated=0;
handles.validation_mode=0;
handles.new_bckg = 1;
guidata(hObject,handles);
update_DA(handles);


% --- Executes during object creation, after setting all properties.
function bckg_poly_order_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bckg_poly_order (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function bckg_poly_order_Callback(hObject, eventdata, handles)
% hObject    handle to bckg_poly_order (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of bckg_poly_order as text
%        str2double(get(hObject,'String')) returns contents of bckg_poly_order as a double
[v,handles]=edit_update(handles,hObject,0,15,3,'%d',1);
if v~=handles.poly_order
	handles.poly_order=v;
    handles.updated=0;
    handles.validation_mode=0;
end
handles.new_bckg = 1;
% Update handles structure
guidata(hObject, handles);
update_DA(handles);


% --- Executes on button press in bckg_save.
function bckg_save_Callback(hObject, eventdata, handles)
% hObject    handle to bckg_save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~isfield(handles,'source_file'), set(handles.status_line,'String','### Load data file ###'); return; end
poly=handles.polynomial;
[fname,pname]=uiputfile('*.dat','Save background polynomial');
fname=strtok(fname,'.');
fullname=[pname fname '.dat'];
save(fullname,'poly','-ascii');

% --- Executes on button press in bckg_load.
function bckg_load_Callback(hObject, eventdata, handles)
% hObject    handle to bckg_load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[fname,pname]=uigetfile('*.dat','Load background polynomial');
fullname=[pname fname];
poly=load(fullname);
handles.polynomial=poly;
handles.updated=0;
handles.validation_mode=0;
handles.new_bckg = 1;
guidata(hObject,handles);
update_DA(handles);


% --- Executes on button press in bckg_add.
function bckg_add_Callback(hObject, eventdata, handles)
% hObject    handle to bckg_add (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

poly0=handles.polynomial;
[fname,pname]=uigetfile('*.dat','Add background polynomial');
fullname=[pname fname];
poly=load(fullname);
weight_str=inputdlg('Weighting','Weighting of new component');
[weight,OK]=str2num(weight_str{1});
if ~OK
    weight=1;
    set(handles.status_line,'String','### Not a number. Defaulting to weighting 1 ###');
end
n0=length(poly0);
n=length(poly);
nn=max([n n0]);
poly1=zeros(1,nn);
for k=1:n0
    poly1(nn+1-k)=poly1(nn+1-k)+poly0(n0+1-k);
end
for k=1:n
    poly1(nn+1-k)=poly1(nn+1-k)+weight*poly(n+1-k);
end
handles.polynomial=poly1;
handles.updated=0;
handles.validation_mode=0;
handles.new_bckg = 1;
guidata(hObject,handles);
update_DA(handles);

% --- Executes on button press in bckg_exp.
function bckg_exp_Callback(hObject, eventdata, handles)
% hObject    handle to bckg_exp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of bckg_exp
flag=get(hObject,'Value');
if flag
    set(handles.bckg_homogeneous,'Value',0);
    set(handles.bckg_none,'Value',0);
    set(handles.bckg_poly,'Value',0);
    set(handles.radiobutton_deernet_bckg,'Value',0);
end
handles.new_bckg = 1;
guidata(hObject,handles);
update_DA(handles);


% --- Executes on button press in fit_Tikhonov.
function fit_Tikhonov_Callback(hObject, eventdata, handles)
% hObject    handle to fit_Tikhonov (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~isfield(handles,'source_file'), set(handles.status_line,'String','### Load data file ###'); return; end
flag=get(handles.select_Tikhonov,'Value');
if flag
	handles.Tikh_nfp_vector=[];
	handles=fit_Tikhonov_new(handles);
	handles.saved=0;
    handles.model_updated=0;
    handles.Tikh_updated=1;
else
    set(handles.status_line,'String','### Select Tikhonov regularization before attempting to fit ###');
end
% Update handles structure
guidata(hObject, handles);
update_DA(handles);


% --- Executes on button press in fit_model.
function fit_model_Callback(hObject, eventdata, handles)
% hObject    handle to fit_model (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~isfield(handles,'source_file'), set(handles.status_line,'String','### Load data file ###'); return; end
flag=get(handles.select_model,'Value');
if flag
    set(handles.error_estimate,'Value',0);
	handles=fit_user_model(handles);
	handles.saved=0;
    handles.model_updated=1;
    handles.Tikh_updated=0;
else
    set(handles.status_line,'String','### Select Model fit before attempting to fit ###');
    handles.model_updated=0;
end
% Update handles structure
guidata(hObject, handles);
update_DA(handles);


% --- Executes during object creation, after setting all properties.
function num_spins_CreateFcn(hObject, eventdata, handles)
% hObject    handle to num_spins (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function num_spins_Callback(hObject, eventdata, handles)
% hObject    handle to num_spins (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of num_spins as text
%        str2double(get(hObject,'String')) returns contents of num_spins as a double

d0=handles.calib_numspins;
v0=handles.n_spins;
lambda=1-exp(handles.calib_numspins*(handles.n_spins-1));
[v,handles]=edit_update(handles,hObject,1,100,handles.n_spins,'%5.2f',0);
handles.n_spins=v;
handles.calib_numspins=log(1-lambda)/(v-1);
set(hObject,'ForegroundColor','g');
guidata(hObject,handles);
update_DA(handles);


% --- Executes on button press in autosave.
function autosave_Callback(hObject, eventdata, handles)
% hObject    handle to autosave (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of autosave


% --- Executes on button press in residual.
function residual_Callback(hObject, eventdata, handles)
% hObject    handle to residual (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of residual
% Update handles structure
guidata(hObject, handles);
update_DA(handles);



% --- Executes on button press in error_estimate.
function error_estimate_Callback(hObject, eventdata, handles)
% hObject    handle to error_estimate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of error_estimate
% Update handles structure
guidata(hObject, handles);
update_DA(handles);




% --- Executes on button press in log_mode.
function log_mode_Callback(hObject, eventdata, handles)
% hObject    handle to log_mode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of log_mode
handles.updated=0;
handles.validation_mode=0;
flag=get(hObject,'Value');
handles.log_Tikh=flag;
guidata(hObject,handles);
update_DA(handles);





% --- Executes on button press in validate_Tikhonov.
function validate_Tikhonov_Callback(hObject, eventdata, handles)
% hObject    handle to validate_Tikhonov (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global main_handles

bmode = get_bckg_mode(handles);

if strcmp(bmode,'d')
    set(handles.main_figure,'Pointer','watch');
    handles = Tikhonov_uncertainty_DEERNet_bckg(handles);
    set(handles.main_figure,'Pointer','arrow');
else
    main_handles=handles;
    % save validation_data main_handles
    h=Tikhonov_validation;
    waitfor(h);
    load('validation_result');
    if ~cancelled
        handles.A_tdip=tdip;
        handles.mean_distr=A_distr;
        handles.A_distr=A_distr;
        handles.A_r=A_r;
        handles.A_sim=A_sim;
        handles.Lcurve_distr=handles.A_distr;
        handles.Lcurve_sim=handles.A_sim;
        handles.bckg_dens=dens;
        handles.bckg_start=bckg_start;
        set(handles.bckg_edit,'String',sprintf('%i',bckg_start));
        handles.A_depth=depth;
        handles.A_dipevo=dipevo;
        handles.A_cluster=clusterp;
        handles.A_bckg=bckg;
        handles.hom_dim=hom_dim;
        handles.moddepth_suppression=moddepth_suppression;
        handles.A_low=dlow;
        handles.A_high=dhigh;
        handles.mask=ones(size(handles.A_distr));
        handles.updated=1;
        handles.validation_mode=1;
        set(handles.error_estimate,'Value',1);
    else
        handles.updated=0;
        handles.validation_mode=0;
        set(handles.error_estimate,'Value',0);
    end
end
guidata(hObject,handles);
update_DA(handles);


% --- Executes on button press in dt_minus.
function dt_minus_Callback(hObject, eventdata, handles)
% hObject    handle to dt_minus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~isfield(handles,'source_file'), set(handles.status_line,'String','### Load data file ###'); return; end
dt=handles.dt-handles.time_grain;
if dt>=handles.min_dt
    handles.dt=dt;
end
guidata(hObject,handles);
update_DA(handles);


% --- Executes on button press in dt_plus.
function dt_plus_Callback(hObject, eventdata, handles)
% hObject    handle to dt_plus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~isfield(handles,'source_file'), set(handles.status_line,'String','### Load data file ###'); return; end
dt=handles.dt+handles.time_grain;
if dt<=handles.max_dt
    handles.dt=dt;
end
guidata(hObject,handles);
update_DA(handles);




% --- Executes on button press in manual_bckg.
function manual_bckg_Callback(hObject, eventdata, handles)
% hObject    handle to manual_bckg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of manual_bckg
handles.updated=0;
if get(hObject,'Value') && isfield(handles,'source_file')
    handles.man_k=handles.density_value/handles.calib_density;
    pstr=sprintf('%6.3f',handles.man_k);
    set(handles.man_bckg_k,'String',pstr);
    handles.man_depth=handles.A_depth;
    pstr=sprintf('%6.3f',handles.A_depth);
    set(handles.man_bckg_depth,'String',pstr);
    set(handles.radiobutton_deernet_bckg,'Value',0);
    handles = set_bckg_mode(handles,handles.old_bckg_mode);
    handles.new_bckg = 1;
end
% Update handles structure
guidata(hObject, handles);
update_DA(handles);



function man_bckg_k_Callback(hObject, eventdata, handles)
% hObject    handle to man_bckg_k (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of man_bckg_k as text
%        str2double(get(hObject,'String')) returns contents of man_bckg_k as a double
[v,handles]=edit_update(handles,hObject,0.0,20,0.2,'%5.3f',0);
if v~=handles.man_depth
	handles.man_k=v;
    handles.updated=0;
    handles.validation_mode=0;
end
handles.new_bckg = 1;
% Update handles structure
guidata(hObject, handles);
update_DA(handles);

% --- Executes during object creation, after setting all properties.
function man_bckg_k_CreateFcn(hObject, eventdata, handles)
% hObject    handle to man_bckg_k (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function man_bckg_depth_Callback(hObject, eventdata, handles)
% hObject    handle to man_bckg_depth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of man_bckg_depth as text
%        str2double(get(hObject,'String')) returns contents of man_bckg_depth as a double
[v,handles]=edit_update(handles,hObject,0.01,0.99,0.33,'%5.3f',0);
if v~=handles.man_depth
    handles.man_depth=v;
    handles.updated=0;
    handles.validation_mode=0;
    handles.new_bckg = 1;
end
% Update handles structure
guidata(hObject, handles);
update_DA(handles);

% --- Executes during object creation, after setting all properties.
function man_bckg_depth_CreateFcn(hObject, eventdata, handles)
% hObject    handle to man_bckg_depth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in man_k_minus.
function man_k_minus_Callback(hObject, eventdata, handles)
% hObject    handle to man_k_minus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[v,handles]=shift_cursor(handles,handles.man_bckg_k,0,20,'%6.3f',handles.man_k,-0.02);
if v~=handles.man_k
    handles.man_k=v;
    handles.updated=0;
    handles.validation_mode=0;
    handles.new_bckg = 1;
end
% Update handles structure
guidata(hObject, handles);
update_DA(handles);


% --- Executes on button press in man_k_plus.
function man_k_plus_Callback(hObject, eventdata, handles)
% hObject    handle to man_k_plus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[v,handles]=shift_cursor(handles,handles.man_bckg_k,0,20,'%6.3f',handles.man_k,0.02);
if v~=handles.man_k
    handles.man_k=v;
    handles.updated=0;
    handles.validation_mode=0;
    handles.new_bckg = 1;
end
% Update handles structure
guidata(hObject, handles);
update_DA(handles);


% --- Executes on button press in man_depth_minus.
function man_depth_minus_Callback(hObject, eventdata, handles)
% hObject    handle to man_depth_minus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[v,handles]=shift_cursor(handles,handles.man_bckg_depth,0.01,0.99,'%6.3f',handles.man_depth,-0.02);
if v~=handles.man_depth
    handles.man_depth=v;
    handles.updated=0;
    handles.validation_mode=0;
    handles.new_bckg = 1;
end
% Update handles structure
guidata(hObject, handles);
update_DA(handles);


% --- Executes on button press in man_depth_plus.
function man_depth_plus_Callback(hObject, eventdata, handles)
% hObject    handle to man_depth_plus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[v,handles]=shift_cursor(handles,handles.man_bckg_depth,0.01,0.99,'%6.3f',handles.man_depth,0.02);
if v~=handles.man_depth
    handles.man_depth=v;
    handles.updated=0;
    handles.validation_mode=0;
    handles.new_bckg = 1;
end
% Update handles structure
guidata(hObject, handles);
update_DA(handles);




% --- Executes on button press in opt_bckg_rms.
function opt_bckg_rms_Callback(hObject, eventdata, handles)
% hObject    handle to opt_bckg_rms (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~isfield(handles,'source_file'), set(handles.status_line,'String','### Load data file ###'); return; end
flag=get(handles.manual_bckg,'Value');
if ~flag
    return;
end
set(handles.main_figure,'Pointer','watch');
drawnow;
[k,depth]=get_manual_bckg(handles);
handles.man_k=k;
pstr=sprintf('%6.3f',k);
set(handles.man_bckg_k,'String',pstr);
handles.man_depth=depth;
pstr=sprintf('%6.3f',depth);
set(handles.man_bckg_depth,'String',pstr);
set(handles.main_figure,'Pointer','arrow');
set(handles.status_line,'String','Ready.');
handles.updated=0;
handles.validation_mode=0;
handles.new_bckg = 1;
% Update handles structure
guidata(hObject, handles);
update_DA(handles);




% --- Executes on button press in renormalize.
function renormalize_Callback(hObject, eventdata, handles)
% hObject    handle to renormalize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of renormalize

guidata(hObject, handles);
update_DA(handles);


% --- Executes on button press in checkbox_guidance.
function checkbox_guidance_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_guidance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_guidance

guidata(hObject, handles);
update_DA(handles);


% --- Executes when user attempts to close main_figure.
function main_figure_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to main_figure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure

shh = get(0,'ShowHiddenHandles');
set(0,'ShowHiddenHandles','on');
delete(findobj(get(0,'Children'),'flat','Tag','StatusBar'));
set(0,'ShowHiddenHandles',shh);
if handles.saved
    if ishandle(handles.original_fig)
        delete(handles.original_fig);
    end
    if ishandle(handles.dipolar_fig)
        delete(handles.dipolar_fig);
    end
    if ishandle(handles.distribution_fig)
        delete(handles.distribution_fig);
    end
	delete(hObject);
else
   answer=questdlg('Do you want to save now?','Fit data not saved');
   if strncmp(answer,'Yes',2)
		return_path=pwd;
		cd(handles.project_dir);
		suggestion=[handles.bas_name '_res.txt']; 
		[fname,pname]=uiputfile('*.txt','Save results',suggestion);
        if isequal(fname,0) || isequal(pname,0)
            return
        end
		% Remove (last) extension, if any
		s=strfind(fname,'.');
		if ~isempty(s)
			fname=fname(1:s(length(s))-1);
    end
		% Remove suffix '_res', if present
		s=strfind(fname,'_res');
		if ~isempty(s)
			fname=fname(1:s(length(s))-1);
    end
		save_result(handles,fname,pname);
		cd(return_path);
    	delete(hObject);
   end
   if strncmp(answer,'No',2)
        if ishandle(handles.original_fig)
            delete(handles.original_fig);
        end
        if ishandle(handles.dipolar_fig)
            delete(handles.dipolar_fig);
        end
        if ishandle(handles.distribution_fig)
            delete(handles.distribution_fig);
        end
       delete(hObject);
   end
end



% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over text_download.
function text_download_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to text_download (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

web('http://www.epr.ethz.ch/software.html','-browser');



% --- Executes on button press in checkbox_model_fit_depth.
function checkbox_model_fit_depth_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_model_fit_depth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_model_fit_depth


% --- Executes on button press in checkbox_ghost.
function checkbox_ghost_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_ghost (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_ghost
% Update handles structure
handles.updated=0;
handles.validation_mode=0;
guidata(hObject, handles);
update_DA(handles);


function edit_oligomerization_Callback(hObject, eventdata, handles)
% hObject    handle to edit_oligomerization (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_oligomerization as text
%        str2double(get(hObject,'String')) returns contents of edit_oligomerization as a double

[v,handles]=edit_update(handles,hObject,2,100,3,'%d',1);
handles.spins_per_object=v;
handles.updated=0;
handles.validation_mode=0;
% Update handles structure
guidata(hObject, handles);
update_DA(handles);

% --- Executes during object creation, after setting all properties.
function edit_oligomerization_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_oligomerization (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on button press in radiobutton_DeerAnalysis.
function radiobutton_DeerAnalysis_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_DeerAnalysis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton_DeerAnalysis

if get(hObject,'Value')
    set(handles.checkbox_no_analysis,'Value',1);
    set(handles.checkbox_no_analysis,'Enable','on');
    set(handles.format_elexsys,'Value',0);
    set(handles.format_winepr,'Value',0);
    set(handles.format_ascii,'Value',0);
else
    set(handles.checkbox_no_analysis,'Enable','off');
end
% Update handles structure
guidata(hObject, handles);
update_DA(handles);

% --- Executes on button press in checkbox_no_analysis.
function checkbox_no_analysis_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_no_analysis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_no_analysis

% Update handles structure
guidata(hObject, handles);
update_DA(handles);


% --- Executes on button press in checkbox_statistics.
function checkbox_statistics_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_statistics (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_statistics


% --- Executes on button press in pushbutton_detach_original.
function pushbutton_detach_original_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_detach_original (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global hMain

if handles.original_attached
    handles.original_attached = false;
    set(hObject,'TooltipString','Attach original data plot');
    handles.original_fig = figure('NumberTitle','off','Name','Original Data','Units','pixels');
    set(hObject,'CData',handles.attach_icon);
    figpos = get(handles.original_fig,'Position');
    pos = get(handles.original_data,'Position');
    handles.original_pos = pos;
    set(handles.original_data,'Parent', handles.original_fig);
    set(handles.original_data,'Units', 'pixels');
    set(handles.original_fig,'CloseRequestFcn',@attach_original);
    set(handles.original_fig,'ResizeFcn',@original_plot_ResizeFcn);
    pos(1)= handles.margin+handles.add_margin;
    pos(2)= handles.margin+handles.add_margin;
    pos(3)= figpos(3)-2*handles.margin-handles.add_margin;
    pos(4)= figpos(4)-2*handles.margin-handles.add_margin;
    set(handles.original_data,'Position',pos);
    hMain = handles;
else
    handles.original_attached = true;
    set(hObject,'TooltipString','Detach original data plot');
    set(hObject,'CData',handles.detach_icon);
    set(handles.original_data,'Parent', handles.main_figure);
    set(handles.original_data,'Units', 'normalized');
    set(handles.original_data,'Position', handles.original_pos);
    delete(handles.original_fig);
end
% Update handles structure
guidata(hObject, handles);
update_DA(handles);

function attach_original(hObject,callbackdata)

global hMain

hMain.original_attached = true;
if ishandle(hMain.pushbutton_detach_original)
    set(hMain.pushbutton_detach_original,'CData',hMain.detach_icon);
    set(hMain.original_data,'Parent', hMain.main_figure);
    set(hMain.original_data,'Units', 'normalized');
    set(hMain.original_data,'Position', hMain.original_pos);
    guidata(hMain.pushbutton_detach_original, hMain);
end
delete(hObject);

% --- Executes on button press in pushbutton_detach_distribution.
function pushbutton_detach_distribution_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_detach_distribution (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global hMain

if handles.distribution_attached
    handles.distribution_attached = false;
    set(hObject,'TooltipString','Attach distance distribution plot');
    handles.distribution_fig = figure('NumberTitle','off','Name','Distance distribution/ L curve','Units','pixels');
    set(hObject,'CData',handles.attach_icon);
    figpos = get(handles.distribution_fig,'Position');
    pos = get(handles.distance_distribution,'Position');
    handles.distribution_pos = pos;
    set(handles.distance_distribution,'Parent', handles.distribution_fig);
    set(handles.distance_distribution,'Units', 'pixels');
    set(handles.distribution_fig,'CloseRequestFcn',@attach_distribution);
    set(handles.distribution_fig,'ResizeFcn',@distribution_plot_ResizeFcn);
    pos(1)= handles.margin+handles.add_margin;
    pos(2)= handles.margin+handles.add_margin;
    pos(3)= figpos(3)-2*handles.margin-handles.add_margin;
    pos(4)= figpos(4)-2*handles.margin-handles.add_margin;
    set(handles.distance_distribution,'Position',pos);
    hMain = handles;
else
    handles.distribution_attached = true;
    set(hObject,'TooltipString','Detach distance distribution plot');
    set(hObject,'CData',handles.detach_icon);
    set(handles.distance_distribution,'Parent', handles.main_figure);
    set(handles.distance_distribution,'Units', 'normalized');
    set(handles.distance_distribution,'Position', handles.distribution_pos);
    delete(handles.distribution_fig);
end
% Update handles structure
guidata(hObject, handles);
update_DA(handles);

function attach_distribution(hObject,callbackdata)

global hMain

hMain.distribution_attached = true;
if ishandle(hMain.pushbutton_detach_distribution)
    set(hMain.pushbutton_detach_distribution,'CData',hMain.detach_icon);
    set(hMain.distance_distribution,'Parent', hMain.main_figure);
    set(hMain.distance_distribution,'Units', 'normalized');
    set(hMain.distance_distribution,'Position', hMain.distribution_pos);
    guidata(hMain.pushbutton_detach_distribution, hMain);
end
delete(hObject);



% --- Executes on button press in pushbutton_detach_dipolar.
function pushbutton_detach_dipolar_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_detach_dipolar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global hMain

if handles.dipolar_attached
    handles.dipolar_attached = false;
    set(hObject,'TooltipString','Attach dipolar evolution plot');
    handles.dipolar_fig = figure('NumberTitle','off','Name','Dipolar evolution function/ spectrum','Units','pixels');
    set(hObject,'CData',handles.attach_icon);
    figpos = get(handles.dipolar_fig,'Position');
    pos = get(handles.dipolar_evolution,'Position');
    handles.dipolar_pos = pos;
    set(handles.dipolar_evolution,'Parent', handles.dipolar_fig);
    set(handles.dipolar_evolution,'Units', 'pixels');
    set(handles.dipolar_fig,'CloseRequestFcn',@attach_dipolar);
    set(handles.dipolar_fig,'ResizeFcn',@dipolar_plot_ResizeFcn);
    pos(1)= handles.margin+handles.add_margin;
    pos(2)= handles.margin+handles.add_margin;
    pos(3)= figpos(3)-2*handles.margin-handles.add_margin;
    pos(4)= figpos(4)-2*handles.margin-handles.add_margin;
    set(handles.dipolar_evolution,'Position',pos);
    hMain = handles;
else
    handles.dipolar_attached = true;
    set(hObject,'TooltipString','Detach dipolar evolution plot');
    set(hObject,'CData',handles.detach_icon);
    set(handles.dipolar_evolution,'Parent', handles.main_figure);
    set(handles.dipolar_evolution,'Units', 'normalized');
    set(handles.dipolar_evolution,'Position', handles.dipolar_pos);
    delete(handles.dipolar_fig);
end
% Update handles structure
guidata(hObject, handles);
update_DA(handles);

function attach_dipolar(hObject,callbackdata)

global hMain

hMain.dipolar_attached = true;
if ishandle(hMain.pushbutton_detach_dipolar)
    set(hMain.pushbutton_detach_dipolar,'CData',hMain.detach_icon);
    set(hMain.dipolar_evolution,'Parent', hMain.main_figure);
    set(hMain.dipolar_evolution,'Units', 'normalized');
    set(hMain.dipolar_evolution,'Position', hMain.dipolar_pos);
    guidata(hMain.pushbutton_detach_dipolar, hMain);
end
delete(hObject);



% --- Executes when model_plot is resized.
function original_plot_ResizeFcn(hObject, callbackdata)
% hObject    handle to model_plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global hMain

newpos=get(hObject,'Position');
newpos(1)= hMain.margin + hMain.add_margin;
newpos(2)= hMain.margin + hMain.add_margin;
newpos(3)= newpos(3)-2*hMain.margin - hMain.add_margin;
newpos(4)= newpos(4)-2*hMain.margin - hMain.add_margin;
set(hMain.original_data,'Position',newpos);

% --- Executes when model_plot is resized.
function distribution_plot_ResizeFcn(hObject, callbackdata)
% hObject    handle to model_plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global hMain

newpos=get(hObject,'Position');
newpos(1)= hMain.margin + hMain.add_margin;
newpos(2)= hMain.margin + hMain.add_margin;
newpos(3)= newpos(3)-2*hMain.margin - hMain.add_margin;
newpos(4)= newpos(4)-2*hMain.margin - hMain.add_margin;
set(hMain.distance_distribution,'Position',newpos);

% --- Executes when model_plot is resized.
function dipolar_plot_ResizeFcn(hObject, callbackdata)
% hObject    handle to model_plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global hMain

newpos=get(hObject,'Position');
newpos(1)= hMain.margin + hMain.add_margin;
newpos(2)= hMain.margin + hMain.add_margin;
newpos(3)= newpos(3)-2*hMain.margin - hMain.add_margin;
newpos(4)= newpos(4)-2*hMain.margin - hMain.add_margin;
set(hMain.dipolar_evolution,'Position',newpos);

function handles = increaseFontSizesIfReq(handles)
% make all fonts smaller on a non-mac-osx computer
persistent fontSizeIncreased
fontSizeIncreased = [];
if ismac()
    % MAC OSX detected; increase font sizes
    if isempty(fontSizeIncreased)
        for afield = fieldnames(handles)'
            afield = afield{1}; %#ok<FXSET>
            try %#ok<TRYNC>
                set(handles.(afield),'FontSize',get(handles.(afield),'FontSize')*(4/3)); % increase font size
            end
        end
        fontSizeIncreased=1; % do not perform this step again.
    end
end

function CData = get_detach_icon

CData(:,:,1) = [   10    0    0    0    0    0    0    0    0    0    0    0    0    0    0   15; ...
   45  220  218  215  213  210  208  204  199  195  191  185  180  173  167   52; ...
   72  255  255  230  195  216  229  244  255  255  255  255  255  255  255   67; ...
   92  255  230   37   31   12    2   95  255  255  255  255  255  255  255   59; ...
  110  255  222    2   83  240  238  241  255  255  255  255  255  255  255   50; ...
  104  255  225    4   61   82  255  255  255  255  255  255  255  255  255   42; ...
   86  255  227    3  237   67   77  238  255  255  255  255  255  255  255   35; ...
   64  255  229    2  255  238   73   71  237  255  255  255  255  255  255   29; ...
   43  255  249  182  255  255  238   80   64  236  255  255  255  255  255   23; ...
   25  255  255  255  255  255  255  255   85   59  235  255  255  255  255   17; ...
   11  255  255  255  255  255  255  255  255   90   52  234  255  255  255   13; ...
    3  255  255  255  255  255  255  255  255  255   96   47  233  255  255    9; ...
    0  255  255  255  255  255  255  255  255  255  255  101   41  244  255    6; ...
    2  255  255  255  255  255  255  255  255  255  255  255  105  154  255    3; ...
    9  255  255  255  255  255  255  255  255  255  255  255  255  189  255   33; ...
   48   34   47   61   76   91  103  113  104   92   79   64   48   35   24   14];

CData(:,:,2) = [   10    0    0    0    0    0    0    0    0    0    0    0    0    0    0   15; ...
   45  220  218  215  213  210  208  204  199  195  191  185  180  173  167   52; ...
   72  255  255  230  195  216  229  244  255  255  255  255  255  255  255   67; ...
   92  255  230   37   31   12    2   95  255  255  255  255  255  255  255   59; ...
  110  255  222    2   83  240  238  241  255  255  255  255  255  255  255   50; ...
  104  255  225    4   61   82  255  255  255  255  255  255  255  255  255   42; ...
   86  255  227    3  237   67   77  238  255  255  255  255  255  255  255   35; ...
   64  255  229    2  255  238   73   71  237  255  255  255  255  255  255   29; ...
   43  255  249  182  255  255  238   80   64  236  255  255  255  255  255   23; ...
   25  255  255  255  255  255  255  255   85   59  235  255  255  255  255   17; ...
   11  255  255  255  255  255  255  255  255   90   52  234  255  255  255   13; ...
    3  255  255  255  255  255  255  255  255  255   96   47  233  255  255    9; ...
    0  255  255  255  255  255  255  255  255  255  255  101   41  244  255    6; ...
    2  255  255  255  255  255  255  255  255  255  255  255  105  154  255    3; ...
    9  255  255  255  255  255  255  255  255  255  255  255  255  189  255   33; ...
   48   34   47   61   76   91  103  113  104   92   79   64   48   35   24   14];

CData(:,:,3) = [   10    0    0    0    0    0    0    0    0    0    0    0    0    0    0   15; ...
   45  220  218  215  213  210  208  204  199  195  191  185  180  173  167   52; ...
   72  255  255  230  195  216  229  244  255  255  255  255  255  255  255   67; ...
   92  255  230   37   31   12    2   95  255  255  255  255  255  255  255   59; ...
  110  255  222    2   83  240  238  241  255  255  255  255  255  255  255   50; ...
  104  255  225    4   61   82  255  255  255  255  255  255  255  255  255   42; ...
   86  255  227    3  237   67   77  238  255  255  255  255  255  255  255   35; ...
   64  255  229    2  255  238   73   71  237  255  255  255  255  255  255   29; ...
   43  255  249  182  255  255  238   80   64  236  255  255  255  255  255   23; ...
   25  255  255  255  255  255  255  255   85   59  235  255  255  255  255   17; ...
   11  255  255  255  255  255  255  255  255   90   52  234  255  255  255   13; ...
    3  255  255  255  255  255  255  255  255  255   96   47  233  255  255    9; ...
    0  255  255  255  255  255  255  255  255  255  255  101   41  244  255    6; ...
    2  255  255  255  255  255  255  255  255  255  255  255  105  154  255    3; ...
    9  255  255  255  255  255  255  255  255  255  255  255  255  189  255   33; ...
   48   34   47   61   76   91  103  113  104   92   79   64   48   35   24   14];

function CData = get_attach_icon

CData(:,:,1) =[    10    0    0    0    0    0    0    0    0    0    0    0    0    0    0   15; ...
   45  220  219  218  216  213  211  207  202  195  191  185  180  173  167   52; ...
   72  255  152  252  243  247  250  253  255  255  255  255  255  255  255   67; ...
   92  255  101  100  255  255  254  254  255  255  255  255  255  255  255   59; ...
  110  255  235   49   67  234  249  251  255  255  255  255  255  255  255   50; ...
  104  255  255  211   66   46  237  255  255  255  255  255  255  255  255   42; ...
   86  255  255  222  226   78   40  219  255  255  255  255  255  255  255   35; ...
   64  255  254  222  255  241   86   32  215  255  255  255  255  255  255   29; ...
   43  255  254  245  255  255  240   90   25  210  255  242   79  157  255   23; ...
   25  255  255  255  255  255  255  255   91   18  203  241   53  172  255   17; ...
   11  255  255  255  255  255  255  255  255   92   11  185   36  190  255   13; ...
    3  255  255  255  255  255  255  255  255  255   95    5   18  205  255    9; ...
    0  255  255  255  255  255  255  255  249  241  241   88    4  211  255    6; ...
    2  255  255  255  255  255  255  255  162   65   80   88    0  179  255    3; ...
    9  255  255  255  255  255  255  255  219  160  146  131  171  196  255   33; ...
   48   34   47   61   76   91  103  113  104   92   79   64   48   35   24   14];

CData(:,:,2) = [   10    0    0    0    0    0    0    0    0    0    0    0    0    0    0   15; ...
   45  220  219  218  216  213  211  207  202  195  191  185  180  173  167   52; ...
   72  255  152  252  243  247  250  253  255  255  255  255  255  255  255   67; ...
   92  255  101  100  255  255  254  254  255  255  255  255  255  255  255   59; ...
  110  255  235   49   67  234  249  251  255  255  255  255  255  255  255   50; ...
  104  255  255  211   66   46  237  255  255  255  255  255  255  255  255   42; ...
   86  255  255  222  226   78   40  219  255  255  255  255  255  255  255   35; ...
   64  255  254  222  255  241   86   32  215  255  255  255  255  255  255   29; ...
   43  255  254  245  255  255  240   90   25  210  255  242   79  157  255   23; ...
   25  255  255  255  255  255  255  255   91   18  203  241   53  172  255   17; ...
   11  255  255  255  255  255  255  255  255   92   11  185   36  190  255   13; ...
    3  255  255  255  255  255  255  255  255  255   95    5   18  205  255    9; ...
    0  255  255  255  255  255  255  255  249  241  241   88    4  211  255    6; ...
    2  255  255  255  255  255  255  255  162   65   80   88    0  179  255    3; ...
    9  255  255  255  255  255  255  255  219  160  146  131  171  196  255   33; ...
   48   34   47   61   76   91  103  113  104   92   79   64   48   35   24   14];

CData(:,:,3) = [   10    0    0    0    0    0    0    0    0    0    0    0    0    0    0   15; ...
   45  220  219  218  216  213  211  207  202  195  191  185  180  173  167   52; ...
   72  255  152  252  243  247  250  253  255  255  255  255  255  255  255   67; ...
   92  255  101  100  255  255  254  254  255  255  255  255  255  255  255   59; ...
  110  255  235   49   67  234  249  251  255  255  255  255  255  255  255   50; ...
  104  255  255  211   66   46  237  255  255  255  255  255  255  255  255   42; ...
   86  255  255  222  226   78   40  219  255  255  255  255  255  255  255   35; ...
   64  255  254  222  255  241   86   32  215  255  255  255  255  255  255   29; ...
   43  255  254  245  255  255  240   90   25  210  255  242   79  157  255   23; ...
   25  255  255  255  255  255  255  255   91   18  203  241   53  172  255   17; ...
   11  255  255  255  255  255  255  255  255   92   11  185   36  190  255   13; ...
    3  255  255  255  255  255  255  255  255  255   95    5   18  205  255    9; ...
    0  255  255  255  255  255  255  255  249  241  241   88    4  211  255    6; ...
    2  255  255  255  255  255  255  255  162   65   80   88    0  179  255    3; ...
    9  255  255  255  255  255  255  255  219  160  146  131  171  196  255   33; ...
   48   34   47   61   76   91  103  113  104   92   79   64   48   35   24   14];


% --- Executes on button press in checkbox_non_negativity.
function checkbox_non_negativity_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_non_negativity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_non_negativity

if get(hObject,'Value')
    handles.non_negativity = true;
else
    handles.non_negativity = false;
end
handles.updated=0;
handles.validation_mode=0;
guidata(hObject, handles);
update_DA(handles);


% --- Executes on mouse press over axes background.
function LCurve_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to distance_distribution (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

pos = get(hObject,'CurrentPoint');
rhoc = pos(1,1);
etac = pos(1,2);
rho = handles.Lcurve_rho;
eta = handles.Lcurve_eta;
rspan = max(rho)-min(rho);
espan = max(eta)-min(eta);
select = 0;
diff = 1e18;
alldiff = zeros(size(eta));
for k = 1:length(rho)
    cdiff = sqrt(((rhoc-rho(k))/rspan)^2 + ((etac-eta(k))/espan)^2);
    alldiff(k) = cdiff;
    if cdiff < diff
        select = k;
        diff = cdiff;
    end
end
handles.regpar = handles.regpar_vector(select); 
handles.regpar_sel = select;
exflag=get(handles.exci_bandwidth_corr,'Value'); % check, if excitation bandwidth correction is selected
set(handles.regpar_edit,'String',num2str(handles.regpar,handles.regpar_edit_strformat));
set(handles.status_line,'String','Recomputing...');
set(handles.main_figure,'Pointer','watch');
drawnow
[r,distr] = get_Tikhonov_new(handles,handles.regpar);
set(handles.status_line,'String','Simulating form factor...');
if exflag
    [sim,sc]=deer_sim(r,distr,handles.A_tdip,handles.bandwidth);
else
    sim=get_td_fit(handles,r,distr);
    sc = 1;
end
handles.moddepth_suppression=sc;
handles.A_sim=sim;
set(handles.status_line,'String','Ready.');
set(handles.main_figure,'Pointer','arrow');
drawnow
handles.A_r=r;
handles.A_distr=distr';
handles.A_low=distr';
handles.A_high=distr';
handles.mask=ones(size(handles.A_distr));
guidata(hObject, handles);
update_DA(handles);

% --- Executes on mouse press over axes background.
function empty_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to distance_distribution (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% --- Executes on button press in checkbox_smooth_scaled.
function checkbox_smooth_scaled_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_smooth_scaled (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_smooth_scaled
guidata(hObject, handles);
update_DA(handles);


% --- Executes on button press in checkbox_deernet.
function checkbox_deernet_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_deernet (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_deernet
guidata(hObject, handles);
update_DA(handles);


% --- Executes on button press in select_deernet.
function select_deernet_Callback(hObject, eventdata, handles)
% hObject    handle to select_deernet (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of select_deernet

flag=get(hObject,'Value');
if flag
    handles.updated=0;
    set(handles.select_Tikhonov,'Value',0);
    set(handles.select_model,'Value',0);
    set(handles.select_APT,'Value',0);
    handles.model_updated=1;
    handles.new_distr = 1;
    handles.validation_mode=0;
    set(handles.L_curve,'Value',0);
    set(handles.L_curve,'Enable','off');
    set(handles.error_estimate,'Value',0);
end
guidata(hObject,handles);
update_DA(handles);


% --- Executes on button press in pushbutton_deernet.
function pushbutton_deernet_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_deernet (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

switch hObject.String
    case 'Compute'
        if ~isfield(handles,'source_file'), set(handles.status_line,'String','### Load data file ###'); return; end
        flag=get(handles.select_deernet,'Value');
        if flag
            handles = compute_deernet(handles);
            handles.saved=0;
            handles.model_updated=0;
            set(handles.error_estimate,'Value',1);
        else
            set(handles.status_line,'String','### Select DEERNet before attempting to fit ###');
        end
    case 'Download'
        button = questdlg('The 2GByte download may take long. Do you want to proceed?','Download neural net sets from server','OK','Cancel','Cancel');
        if strcmp(button,'OK')
            set(handles.main_figure,'Pointer','watch');
            netset_defs = which('distributed_netsets.txt');
            deernet_dir = fileparts(netset_defs);
            existing_nets = 0;
            for nd = 1:length(handles.net_sets)
                set(handles.status_line,'String',sprintf('Downloading %s...',handles.net_sets(nd).directory));
                drawnow
                unzip(handles.net_sets(nd).url,deernet_dir);
                if exist(sprintf('%s%c1.mat',handles.net_sets(nd).directory,filesep),'file')
                    existing_nets = existing_nets + 1;
                    menu{existing_nets} = handles.net_sets(nd).name;
                end
            end
            % menu{existing_nets + 1} = '<check for updates>'; % needs
            % adapatation of path
            set(handles.popupmenu_deernet,'String',menu);
            if existing_nets > 0
                set(handles.pushbutton_deernet,'String','Compute');
                set(handles.pushbutton_deernet,'TooltipString','Starts neural network analysis');
            else
                set(handles.pushbutton_deernet,'String','Download');
                set(handles.pushbutton_deernet,'TooltipString','Downloads neural network sets from server');
            end
            set(handles.status_line,'String','Neural network download completed.');
            set(handles.main_figure,'Pointer','arrow');
            drawnow
        end
end
% Update handles structure
guidata(hObject, handles);
update_DA(handles);

% --- Executes on selection change in popupmenu_deernet.
function popupmenu_deernet_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_deernet (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_deernet contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_deernet

old_net_set = handles.net_set;
choice = get(hObject,'Value');
if length(handles.net_sets) >= choice
    handles.net_set = handles.net_sets(choice).directory;
end
if handles.select_deernet.Value && ~strcmp(old_net_set,handles.net_set)
    handles.updated = 0;
    handles.model_updated = 1;
end
% Update handles structure
guidata(hObject, handles);
update_DA(handles);

% --- Executes during object creation, after setting all properties.
function popupmenu_deernet_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_deernet (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_comparative.
function checkbox_comparative_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_comparative (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_comparative

% Update handles structure
guidata(hObject, handles);
update_DA(handles);


% --- Executes on button press in radiobutton_deernet_bckg.
function radiobutton_deernet_bckg_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_deernet_bckg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton_deernet_bckg
flag=get(hObject,'Value');
if flag
    set_bckg_mode(handles,'d');
end
handles.updated=0;
handles.model_updated=0;
handles.validation_mode=0;
handles.new_bckg = 1;


% Update handles structure
guidata(hObject, handles);
update_DA(handles);

% --- Executes on button press in checkbox_deernet_error_bckg.
function checkbox_deernet_error_bckg_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_deernet_error_bckg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_deernet_error_bckg
% Update handles structure
guidata(hObject, handles);
update_DA(handles);

% --- Executes on button press in checkbox_comparative_bckg.
function checkbox_comparative_bckg_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_comparative_bckg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

guidata(hObject,handles);
update_DA(handles);



% --- Executes on button press in pushbutton_validate_model.
function pushbutton_validate_model_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_validate_model (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

bmode = get_bckg_mode(handles);

if strcmp(bmode,'d')
    set(handles.main_figure,'Pointer','watch');
    handles = user_model_uncertainty_DEERNet_bckg(handles);
    set(handles.main_figure,'Pointer','arrow');
else
    set(hObject,'Enable','off');
end

guidata(hObject,handles);
update_DA(handles);


% --- Executes on button press in regpar_GCV.
function regpar_GCV_Callback(hObject, eventdata, handles)
% hObject    handle to regpar_GCV (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~isfield(handles,'source_file'), set(handles.status_line,'String','### Load data file ###'); return; end
if isempty(handles.Lcurve_sim), set(handles.status_line,'String','### Tikhonov fit required. ###'); return; end
if numel(handles.A_r)~=numel(handles.A_distr), set(handles.status_line,'String','### Tikhonov L curve fit required. ###'); return; end
if numel(handles.regpar_vector)==1, set(handles.status_line,'String','### Tikhonov L curve fit required. ###'); return; end

idxGCV = handles.regpar_opt_GCV;
handles.regpar_sel = idxGCV;
handles.regpar = handles.regpar_vector(idxGCV);

set(handles.regpar_edit,'String',num2str(handles.regpar,handles.regpar_edit_strformat));
set(handles.status_line,'String','Recomputing...');
set(handles.main_figure,'Pointer','watch');
drawnow

[r,distr] = get_Tikhonov_new(handles,handles.regpar);

set(handles.status_line,'String','Simulating form factor...');
% check if excitation bandwidth correction is selected
exBWcorr = get(handles.exci_bandwidth_corr,'Value'); 
if exBWcorr
    [sim,scale] = deer_sim(r,distr,handles.A_tdip,handles.bandwidth);
elseif length(handles.A_tdip) > 1024
    sim = deer_sim(r,distr,handles.A_tdip);
    scale = 1;
else
    sim = get_td_fit(handles,r,distr);
    scale = 1;
end
handles.moddepth_suppression = scale;
handles.A_sim = sim;

set(handles.status_line,'String','Ready.');
set(handles.main_figure,'Pointer','arrow');
drawnow
handles.A_r=r;
handles.A_distr=distr';
handles.A_low=distr';
handles.A_high=distr';
handles.mask=ones(size(handles.A_distr));
% Update handles structure
guidata(hObject, handles);
update_DA(handles);
