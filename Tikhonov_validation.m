function varargout = Tikhonov_validation(varargin)
% TIKHONOV_VALIDATION M-file for Tikhonov_validation.fig
%      TIKHONOV_VALIDATION, by itself, creates a new TIKHONOV_VALIDATION or raises the existing
%      singleton*.
%
%      H = TIKHONOV_VALIDATION returns the handle to a new TIKHONOV_VALIDATION or the handle to
%      the existing singleton*.
%
%      TIKHONOV_VALIDATION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TIKHONOV_VALIDATION.M with the given input arguments.
%
%      TIKHONOV_VALIDATION('Property','Value',...) creates a new TIKHONOV_VALIDATION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Tikhonov_validation_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Tikhonov_validation_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Tikhonov_validation

% Last Modified by GUIDE v2.5 07-Jan-2011 08:08:03

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Tikhonov_validation_OpeningFcn, ...
                   'gui_OutputFcn',  @Tikhonov_validation_OutputFcn, ...
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


% --- Executes just before Tikhonov_validation is made visible.
function Tikhonov_validation_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Tikhonov_validation (see VARARGIN)

global main_handles

% Choose default command line output for Tikhonov_validation
handles.output = hObject;

% load validation_data
handles.calib_density=main_handles.calib_density;
handles.dens=main_handles.bckg_dens;
handles.mod_depth=main_handles.A_depth;
handles.hom_dim=main_handles.hom_dim;
defaults=load('validation_defaults.dat');
if defaults(19)>main_handles.cutoff
    defaults(9)=main_handles.cutoff;
end
handles.ghost_suppression = get(main_handles.checkbox_ghost,'Value');
handles.spins_per_object = main_handles.spins_per_object;
handles.constraint_level=defaults(13);
handles.pruning=defaults(12);
set(handles.prune_level,'String',sprintf('%4.2f',handles.pruning));
handles.factor_noise=defaults(2);
handles.min_cons=1.0;
handles.max_cons=10.0;
handles.upper_c_bound=defaults(14);
set(handles.upper_bound,'String',sprintf('%5.2f',handles.upper_c_bound));
handles.bas_name=main_handles.bas_name;
set(handles.noise_factor,'String',sprintf('%4.1f',handles.factor_noise));
handles.trials_noise=defaults(1);
set(handles.noise_trials,'String',num2str(handles.trials_noise));
handles.trials_bckg=defaults(3);
set(handles.bckg_trials,'String',num2str(handles.trials_bckg));
handles.trials_depth=defaults(6);
set(handles.mod_depth_trials,'String',num2str(handles.trials_depth));
handles.trials_dim=defaults(9);
set(handles.dim_trials,'String',num2str(handles.trials_dim));
handles.trials_bckg_start=defaults(17);
set(handles.edit_bckg_start_trials,'String',num2str(handles.trials_bckg_start));
handles.trials_total=handles.trials_noise*handles.trials_bckg*handles.trials_depth*handles.trials_dim*handles.trials_bckg_start;
set(handles.total_trials,'String',num2str(handles.trials_total));
handles.interactive=defaults(16);
handles.stop_flag=0;
handles.copy_flag=0;
handles.computed=0;
set(handles.error_estimate,'Value',1);
set(handles.prune,'Enable','off');
handles.dt=main_handles.dt;
handles.fit_rms_value=main_handles.fit_rms_value;
pstr=sprintf('%8.6f',handles.fit_rms_value);
set(handles.rms_level,'String',pstr);
handles.min_bckg=defaults(4)/handles.calib_density;
pstr=sprintf('%6.4f',handles.min_bckg*handles.calib_density);
set(handles.bckg_min,'String',pstr);
handles.max_bckg=defaults(5)/handles.calib_density;
pstr=sprintf('%6.4f',handles.max_bckg*handles.calib_density);
set(handles.bckg_max,'String',pstr);
handles.min_dim=defaults(10);
pstr=sprintf('%5.2f',handles.min_dim);
set(handles.dim_min,'String',pstr);
handles.max_dim=defaults(11);
pstr=sprintf('%5.2f',handles.max_dim);
set(handles.dim_max,'String',pstr);
handles.mod_depth_min=defaults(7);
pstr=sprintf('%5.3f',handles.mod_depth_min);
set(handles.min_mod_depth,'String',pstr);
handles.mod_depth_max=defaults(8);
pstr=sprintf('%5.3f',handles.mod_depth_max);
set(handles.max_mod_depth,'String',pstr);
handles.auto_save=defaults(15);
dt = main_handles.texp(2)-main_handles.texp(1);
defaults(18) = 1000*dt*round(0.2*max(main_handles.texp)/dt); % adapts to data length
handles.min_bckg_start=defaults(18);
pstr=sprintf('%i',handles.min_bckg_start);
set(handles.edit_bckg_start_min,'String',pstr);
defaults(19) = 1000*dt*round(0.6*max(main_handles.texp)/dt); % adapts to data length
handles.max_bckg_start=defaults(19);
pstr=sprintf('%i',handles.max_bckg_start);
set(handles.edit_bckg_start_max,'String',pstr);
distr0=main_handles.A_distr;
handles.distr_std=0*distr0;
set(handles.status_line,'String',sprintf('%s%6.3f%s%6.3f','Initial density: ',main_handles.bckg_dens*handles.calib_density,', Initial depth: ',main_handles.A_depth));

handles.theor=main_handles.theor;
if handles.theor
    handles.th_distr=main_handles.th_distr;
    handles.th_r=main_handles.th_r;
end

figname=sprintf('%s%s','Tikhonov validation with regularization parameter ',num2str(main_handles.regpar)); % tell user, which file is current
set(handles.Tikhonov_validation,'Name',figname);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Tikhonov_validation wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Tikhonov_validation_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in noise_test.
function noise_test_Callback(hObject, eventdata, handles)
% hObject    handle to noise_test (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of noise_test
update_trials(handles);



function noise_factor_Callback(hObject, eventdata, handles)
% hObject    handle to noise_factor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of noise_factor as text
%        str2double(get(hObject,'String')) returns contents of noise_factor as a double
[v,handles]=edit_update(handles,hObject,1,10,2,'%5.2f',0);
if v~=handles.factor_noise
	handles.factor_noise=v;
end
% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function noise_factor_CreateFcn(hObject, eventdata, handles)
% hObject    handle to noise_factor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function noise_trials_Callback(hObject, eventdata, handles)
% hObject    handle to noise_trials (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of noise_trials as text
%        str2double(get(hObject,'String')) returns contents of noise_trials as a double
[v,handles]=edit_update(handles,hObject,1,50,5,'%d',1);
if v~=handles.trials_noise
	handles.trials_noise=v;
end;
% Update handles structure
guidata(hObject, handles);
update_trials(handles);


% --- Executes during object creation, after setting all properties.
function noise_trials_CreateFcn(hObject, eventdata, handles)
% hObject    handle to noise_trials (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in background_density.
function background_density_Callback(hObject, eventdata, handles)
% hObject    handle to background_density (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of background_density
update_trials(handles);



function bckg_min_Callback(hObject, eventdata, handles)
% hObject    handle to bckg_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of bckg_min as text
%        str2double(get(hObject,'String')) returns contents of bckg_min as a double
[v,handles]=edit_update(handles,hObject,0,handles.max_bckg*handles.calib_density,0.1,'%6.4f',0);
v=v/handles.calib_density;
if v~=handles.min_bckg
	handles.min_bckg=v;
end
% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function bckg_min_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bckg_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function bckg_max_Callback(hObject, eventdata, handles)
% hObject    handle to bckg_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of bckg_max as text
%        str2double(get(hObject,'String')) returns contents of bckg_max as a double
[v,handles]=edit_update(handles,hObject,handles.min_bckg*handles.calib_density,10,0.2,'%6.4f',0);
v=v/handles.calib_density;
if v~=handles.max_bckg
	handles.max_bckg=v;
end
% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function bckg_max_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bckg_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function bckg_trials_Callback(hObject, eventdata, handles)
% hObject    handle to bckg_trials (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of bckg_trials as text
%        str2double(get(hObject,'String')) returns contents of bckg_trials as a double
[v,handles]=edit_update(handles,hObject,1,50,5,'%d',1);
if v~=handles.trials_bckg
	handles.trials_bckg=v;
end
% Update handles structure
guidata(hObject, handles);
update_trials(handles);


% --- Executes during object creation, after setting all properties.
function bckg_trials_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bckg_trials (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in background_dim.
function background_dim_Callback(hObject, eventdata, handles)
% hObject    handle to background_dim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of background_dim
update_trials(handles);



function dim_trials_Callback(hObject, eventdata, handles)
% hObject    handle to dim_trials (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dim_trials as text
%        str2double(get(hObject,'String')) returns contents of dim_trials as a double
[v,handles]=edit_update(handles,hObject,1,50,5,'%d',1);
if v~=handles.trials_dim
	handles.trials_dim=v;
end
% Update handles structure
guidata(hObject, handles);
update_trials(handles);

% --- Executes during object creation, after setting all properties.
function dim_trials_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dim_trials (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function dim_min_Callback(hObject, eventdata, handles)
% hObject    handle to dim_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dim_min as text
%        str2double(get(hObject,'String')) returns contents of dim_min as a double
[v,handles]=edit_update(handles,hObject,1,handles.max_dim,3,'%5.2f',0);
if v~=handles.min_dim
	handles.min_dim=v;
end
% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function dim_min_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dim_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function dim_max_Callback(hObject, eventdata, handles)
% hObject    handle to dim_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dim_max as text
%        str2double(get(hObject,'String')) returns contents of dim_max as a double
[v,handles]=edit_update(handles,hObject,handles.min_dim,4,3,'%5.2f',0);
if v~=handles.max_dim
	handles.max_dim=v;
end
% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function dim_max_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dim_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in compute.
function compute_Callback(hObject, eventdata, handles)
% hObject    handle to compute (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global main_handles

% load validation_data
ntrials=handles.trials_noise;
if ~get(handles.noise_test,'Value')
    ntrials=1;
    handles.factor_noise=1;
    msg=sprintf('%5.2f',handles.factor_noise);
    set(handles.noise_factor,'String',msg);
end
nlev_vec=handles.factor_noise*ones(1,ntrials);
btrials=handles.trials_bckg;
if ~get(handles.background_density,'Value')
    dens_vec=main_handles.bckg_dens;
    btrials=1;
else
    dens_vec=linspace(handles.min_bckg,handles.max_bckg,btrials);
end
ltrials=handles.trials_depth;
if ~get(handles.vary_depth,'Value')
    depth_vec=main_handles.A_depth;
    ltrials=1;
else
    depth_vec=linspace(handles.mod_depth_min,handles.mod_depth_max,ltrials);
end
dtrials=handles.trials_dim;
if ~get(handles.background_dim,'Value')
    dim_vec=handles.hom_dim;
    dtrials=1;
else
    dim_vec=linspace(handles.min_dim,handles.max_dim,dtrials);
end
strials=handles.trials_bckg_start;
if ~get(handles.checkbox_bckg_start,'Value')
    bckg_start_vec=main_handles.bckg_start;
    strials=1;
else
    bckg_start_vec=linspace(handles.min_bckg_start,handles.max_bckg_start,strials);
end
t_cutoff=main_handles.cutoff;
texp=main_handles.texp;
handles.texp=texp;
vexp=main_handles.vexp;
handles.vexp=vexp;
[~,ztpoi]=min(abs(texp));
ttemp=texp-t_cutoff*ones(size(texp))/1000;
[~,pcutoff]=min(abs(ttemp));
tdip=texp(ztpoi:pcutoff);
handles.A_tdip=tdip;
kc=1;
total_trials=ntrials*btrials*ltrials*dtrials*strials;
if total_trials<5
    set(handles.status_line,'String','### Validation requires at least 5 trial computations. ###');
    return
else
    store_bckg_start=main_handles.bckg_start;
    store_hom_dim=main_handles.hom_dim;
    if ~get(handles.statusbar_off,'Value')
        status_figure('Validation. Close to stop.');
    else
        set(handles.Tikhonov_validation,'Pointer','watch');
        drawnow;
    end;
    poi=0;
    npoi=0;
%     distr0=main_handles.A_distr;
    set(main_handles.select_L_curve,'Value',0);
    trial_bckg=zeros(total_trials,length(tdip));
    trial_dipevo=zeros(total_trials,length(tdip));
    trial_cluster=zeros(total_trials,length(tdip));
    trial_distr=zeros(total_trials,1024);
    trial_sim=zeros(total_trials,length(tdip));
    trial_sim0=zeros(total_trials,length(tdip));
    rms_vec=zeros(1,total_trials);
    dens_vals=zeros(1,total_trials);
    depth_vals=zeros(1,total_trials);
    mod_depths=zeros(1,total_trials);
    dim_vals=zeros(1,total_trials);
    sc_vals=zeros(1,total_trials);
    start_vals=zeros(1,total_trials);
    tic;
    for kn=1:ntrials
        nfac=sqrt(nlev_vec(kn)^2-1);
        noise=nfac*handles.fit_rms_value*randn(size(vexp));
        for kd=1:dtrials
            hom_dim=dim_vec(kd);
            for kb=1:btrials
                dens=dens_vec(kb);
                for kl=1:ltrials
                    mod_depth=depth_vec(kl);
                    for ks=1:strials
                        t_bckg=bckg_start_vec(ks);
                        npoi=npoi+1;
                        if get(handles.checkbox_bckg_start,'Value')
                            t_cutoff=main_handles.cutoff;
                            ttemp=1000*texp-t_bckg*ones(size(texp));
                            [~,nofitp0]=min(abs(ttemp));
                            ttemp=1000*texp-t_cutoff*ones(size(texp));
                            [~,pcutoff]=min(abs(ttemp));
                            main_handles.bckg_start=t_bckg;
                            main_handles.hom_dim=hom_dim;
                            t_fit=texp(nofitp0:pcutoff); % time window of baseline region
                            td_fit=real(vexp(nofitp0:pcutoff)); % experimental data in this window
                            [bckg,main_handles]=fit_bckg(main_handles,texp,t_fit,td_fit);
                            mod_depth = 1-bckg(1);                         
                        else
                            % Background simulation
                            v1=[dens (1-mod_depth)];
                            bckg=decaynD(v1,texp,hom_dim);
                        end
                        % bckg=(1-main_handles.man_depth)*bckg;
                        if handles.interactive
                            axes(handles.axes1);
                            cla;
                            plot(texp,real(vexp),'k');
                            hold on;
                            plot(texp,real(bckg),'r');
                            xlabel('t (탎)');
                            axis([min(texp),max(texp),min([min(bckg) min(real(vexp))]),1.1*max([max(bckg) 1.1*max(real(vexp))])]);
                        end
                        dipevo=real(vexp)+noise-bckg;
                        cluster=real(vexp)+noise;
                        dipevo=dipevo./bckg; % divide by background, eqn [13]
                        cluster=cluster./bckg;
                        if handles.ghost_suppression
                            cluster=cluster.^(1/(handles.spins_per_object-1));
                            dipevo=cluster-ones(size(cluster));
                        end                        
                        cluster=cluster/max(cluster);
                        dipevo=dipevo(ztpoi:pcutoff);
                        cluster=cluster(ztpoi:pcutoff);
                        bckg=bckg(ztpoi:pcutoff);
                        tdip=texp(ztpoi:pcutoff);
                        dipevo=dipevo/max(dipevo);
                        cluster=cluster/max(cluster);
                        main_handles.A_tdip=tdip;
                        main_handles.A_dipevo=dipevo;
                        main_handles.A_cluster=cluster;
                        msg=sprintf('%s%i%s%i%s%5.2f%s%6.4f%s%5.3f%s%5.2f','Trial ',npoi,' of ',total_trials,', Noise: ',nlev_vec(kn),', Density: ',dens_vec(kb)*handles.calib_density,', Mod. depth: ',depth_vec(kl),', Dim.: ',dim_vec(kd));
                        set(handles.param_line,'String',msg);
                        [A_r,distr,sim,rms,sc]=get_Tikhonov(handles,main_handles);
                        if ~isempty(distr)
                            poi=poi+1;
                            modsim=ones(size(sim))-sim;
                            sim1=ones(size(modsim))-sc*modsim;
                            if handles.interactive
                                axes(handles.axes2);
                                cla;
                                plot(tdip,cluster,'k');
                                hold on;
                                plot(tdip,sim1,'r');
                                axis([min(tdip),max(tdip),min([min(sim1) min(cluster)]),max([max(sim1) max(cluster)])]);
                                xlabel('t (탎)');
                                axes(handles.axes3);
                                cla;
                                plot(A_r,distr,'k');
                                axis([min(A_r),max(A_r),-0.1*max(distr),1.1*max(distr)]);
                                xlabel('r (nm)');
                                drawnow
                            end
                            trial_bckg(poi,:)=bckg;
                            trial_dipevo(poi,:)=dipevo;
                            trial_cluster(poi,:)=cluster;
                            trial_distr(poi,1:length(distr))=distr;
                            trial_sim(poi,:)=sim1;
                            trial_sim0(poi,:)=sim;
                            rms_vec(poi)=rms;
                            dens_vals(poi)=dens_vec(kb);
                            depth_vals(poi)=depth_vec(kl);
                            mod_depths(poi) = mod_depth;
                            dim_vals(poi)=dim_vec(kd);
                            start_vals(poi)=bckg_start_vec(ks);
                            sc_vals(poi)=sc;
                        end
                        if ~get(handles.statusbar_off,'Value')
                            comp_status=status_figure(npoi/total_trials);
                            drawnow;
                            if ~comp_status, stop_flag=1; else stop_flag=0; end
                        else
                            stop_flag=0;
                        end
                        if stop_flag==1
                            break;
                        end
                    end
                end
                if stop_flag==1
                    break;
                end
            end
            if stop_flag==1
                break;
            end
        end
        if stop_flag==1
            break;
        end
    end
    if ~get(handles.statusbar_off,'Value')
        if comp_status, status_figure(1); end
    else
        set(handles.Tikhonov_validation,'Pointer','arrow');
        drawnow;
    end
    handles.successful_trials=poi;
    trial_bckg=trial_bckg(1:poi,:);
    trial_dipevo=trial_dipevo(1:poi,:);
    trial_cluster=trial_cluster(1:poi,:);
    trial_distr=trial_distr(1:poi,1:length(distr));
    trial_sim=trial_sim(1:poi,:);
    trial_sim0=trial_sim0(1:poi,:);
    rms_vec=rms_vec(1:poi);
    handles.dens_vals=dens_vals(1:poi);
    handles.depth_vals=depth_vals(1:poi);
    handles.dim_vals=dim_vals(1:poi);
    handles.sc_vals=sc_vals(1:poi);
    handles.start_vals=start_vals(1:poi);
    endtime=toc;
    msg=sprintf('%s%8.1f%s','Total computation time: ',endtime,' s');
    if stop_flag==1
        msg=sprintf('%s%s%i%s',msg,'. Computation stopped after ',npoi,' trials.');
    end
    set(handles.status_line,'String',msg);
    [min_rms,rpoi]=min(rms_vec);
%     tz=tz_vals(rpoi);
%     tz=handles.dt*round(tz/handles.dt);
%     msg=sprintf('%d',tz);
%     set(main_handles.bckg_edit,'String',msg);
    handles.dens=dens_vals(rpoi);
    handles.mod_depth=depth_vals(rpoi);
    handles.bckg_start=start_vals(rpoi);
    msg=sprintf('%5.2f',dim_vals(rpoi));
    set(main_handles.bckg_dim_edit,'String',msg);
    set(main_handles.bckg_fit_dim,'Value',0);
    distr_std=std(trial_distr);
    mean_distr=mean(trial_distr);
    handles.mean_distr=mean_distr;
    handles.distr_std=distr_std;
    mean_depth=mean(mod_depths);
    std_depth=std(mod_depths);
    handles.mean_depth=mean_depth;
    handles.std_depth=std_depth;    
    handles.texp=texp;
    handles.vexp=vexp;
    handles.A_selected=rpoi;
    handles.A_r=A_r;
    handles.A_tdip=tdip;
    handles.A_distr=trial_distr(rpoi,:);
    handles.A_sim=trial_sim0(rpoi,:);
    handles.A_cluster=trial_cluster(rpoi,:);
    handles.A_dipevo=trial_dipevo(rpoi,:);
    handles.A_bckg=trial_bckg(rpoi,:);
    handles.trial_distr=trial_distr;
    handles.trial_sim=trial_sim;
    handles.trial_sim0=trial_sim0;
    handles.trial_cluster=trial_cluster;
    handles.trial_dipevo=trial_dipevo;
    handles.trial_bckg=trial_bckg;
    handles.trial_rmsd=rms_vec;
    % handles.bckg_start=tz;
    handles.hom_dim=dim_vals(rpoi);
    handles.moddepth_suppression=sc_vals(rpoi);
    msg=sprintf('%s%8.6f%s%6.4f%s%5.3f%s%5.2f','Best r.m.s. ',min_rms,' at density ',dens_vals(rpoi)*handles.calib_density,' and mod. depth ',depth_vals(rpoi),' with bckg. dim ',dim_vals(rpoi)); 
    set(handles.param_line,'String',msg);
    set(handles.parameter_set,'String',sprintf('%i',rpoi));
    set(handles.best_rmsd,'String',sprintf('%8.6f',min_rms));
    main_handles.bckg_start=store_bckg_start;
    main_handles.hom_dim=store_hom_dim;
end
if handles.auto_save
    fname=[handles.bas_name '_validation_sets'];
    save(fname,'trial_bckg','trial_dipevo','trial_distr','trial_sim','dens_vals','depth_vals','dim_vals','sc_vals','mean_depth','std_depth');
end
handles.computed=1;
set(handles.prune,'Enable','on');
% Update handles structure
guidata(hObject, handles);
update_plots(handles);

% --- Executes on button press in cancel.
function cancel_Callback(hObject, eventdata, handles)
% hObject    handle to cancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
distr_std=0*handles.distr_std;
cancelled=1;
save validation_result distr_std cancelled
close(Tikhonov_validation);

% --- Executes on button press in close_validation.
function close_validation_Callback(hObject, eventdata, handles)
% hObject    handle to close_validation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~handles.computed % if no result exists, Cancel instead of closing politely
    distr_std=0*handles.distr_std;
    cancelled=1;
    save validation_result distr_std cancelled
    close(Tikhonov_validation);
    return
end
distr_std=handles.distr_std;
cancelled=0;
sel=handles.A_selected;
A_distr = mean(handles.trial_distr);
A_distr(A_distr<0) = 0;
dlow = min(handles.trial_distr);
dhigh = max(handles.trial_distr);
dlow(dlow<0) = 0;
A_sim=handles.trial_sim0(sel,:);
A_r=handles.A_r;
tdip=handles.A_tdip;
bckg=handles.trial_bckg(sel,:);
dipevo=handles.trial_dipevo(sel,:);
clusterp=handles.trial_cluster(sel,:);
dens=handles.dens_vals(sel);
depth=handles.depth_vals(sel);
hom_dim=handles.dim_vals(sel);
bckg_start=handles.start_vals(sel);
moddepth_suppression=handles.moddepth_suppression;

wfile=fopen([handles.bas_name '_validation.txt'],'wt');
fprintf(wfile,'%s%s\n\n','Tikhonov validation for file ',handles.bas_name);
fprintf(wfile,'%s\n\n','Selected data set:');
fprintf(wfile,'%s%6.3f%s%6.3f\n','Density calibrated:',dens*handles.calib_density,' uncalibrated:',dens);
fprintf(wfile,'%s%i ns\n','Background start:',bckg_start);
fprintf(wfile,'%s%6.3f\n','Mod. depth:',depth);
fprintf(wfile,'%s%5.2f\n','Dimension:',hom_dim);
fprintf(wfile,'%s%9.6f\n','r.m.s.d.:',handles.trial_rmsd(sel));
fprintf(wfile,'%s%5.2f%s%5.2f%s\n','Constraints: (',handles.min_cons,', ',handles.max_cons,') nm');
dens_min=min(handles.dens_vals);
dens_max=max(handles.dens_vals);
depth_min=min(handles.depth_vals);
depth_max=max(handles.depth_vals);
dim_min=min(handles.dim_vals);
dim_max=max(handles.dim_vals);
fprintf(wfile,'\n%s\n\n','All data sets:');
fprintf(wfile,'%s%6.3f%s%6.3f%s\n','Background start range: (',min(handles.start_vals),',',max(handles.start_vals),') ns');
fprintf(wfile,'%s%6.3f%s%6.3f%s\n','Density range: (',dens_min*handles.calib_density,',',dens_max*handles.calib_density,')');
fprintf(wfile,'%s%6.3f%s%6.3f%s\n','Mod. depth range: (',depth_min,',',depth_max,')');
fprintf(wfile,'%s%5.2f%s%5.2f%s\n','Dimension range: (',dim_min,',',dim_max,')');
fprintf(wfile,'%s%9.6f','Best r.m.s.d.:',min(handles.trial_rmsd));
fclose(wfile);

save validation_result distr_std cancelled A_r A_distr A_sim bckg dipevo clusterp tdip bckg_start dens depth hom_dim moddepth_suppression dlow dhigh 
close(Tikhonov_validation);

function update_trials(handles)
% Updates the total number of trials
total=1;
if get(handles.noise_test,'Value')
    total=total*handles.trials_noise;
end;
if get(handles.background_density,'Value')
    total=total*handles.trials_bckg;
end;
if get(handles.vary_depth,'Value')
    total=total*handles.trials_depth;
end;
if get(handles.background_dim,'Value')
    total=total*handles.trials_dim;
end;
if get(handles.checkbox_bckg_start,'Value')
    total=total*handles.trials_bckg_start;
end;
handles.trials_total=total;
pstr=sprintf('%d',total);
set(handles.total_trials,'String',pstr);

function [r,distr,sim,rms,sc]=get_Tikhonov(handles,main_handles)
%
% Tikhonov regularization, with or without excitation bandwidth correction
% and with or without computation of a whole L curve
% Matlab-internal version
%
% G. Jeschke, 2015
%


set(handles.Tikhonov_validation,'Pointer','watch');
drawnow;
[r,distr,rho,eta,reg_param,corner] = get_Tikhonov_new(main_handles,main_handles.regpar);
set(handles.Tikhonov_validation,'Pointer','arrow');
drawnow;
exflag=get(main_handles.exci_bandwidth_corr,'Value');
if exflag,
    sim = deer_sim(r,distr,main_handles.A_tdip,main_handles.bandwidth);
else
    sim = get_td_fit(main_handles,r,distr);
end;
sim = real(sim);
modsim=ones(size(sim))-sim;
modexp=ones(size(main_handles.A_cluster))-main_handles.A_cluster;
sc=sum(modexp.*modexp)/sum(modsim.*modexp);
sc = real(sc);
sim1=ones(size(modsim))-sc*modsim;
diff0=main_handles.A_cluster-sim1;
rms=sqrt(sum(diff0.*diff0)/(length(diff0)-1));


% --- Executes on button press in vary_depth.
function vary_depth_Callback(hObject, eventdata, handles)
% hObject    handle to vary_depth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of vary_depth



function min_mod_depth_Callback(hObject, eventdata, handles)
% hObject    handle to min_mod_depth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of min_mod_depth as text
%        str2double(get(hObject,'String')) returns contents of min_mod_depth as a double
[v,handles]=edit_update(handles,hObject,0.01,handles.mod_depth_max,0.2,'%5.3f',0);
if v~=handles.mod_depth_min,
	handles.mod_depth_min=v;
end;
% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function min_mod_depth_CreateFcn(hObject, eventdata, handles)
% hObject    handle to min_mod_depth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function max_mod_depth_Callback(hObject, eventdata, handles)
% hObject    handle to max_mod_depth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of max_mod_depth as text
%        str2double(get(hObject,'String')) returns contents of max_mod_depth as a double
[v,handles]=edit_update(handles,hObject,handles.mod_depth_min,0.9,0.5,'%5.3f',0);
if v~=handles.mod_depth_max,
	handles.mod_depth_max=v;
end;
% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function max_mod_depth_CreateFcn(hObject, eventdata, handles)
% hObject    handle to max_mod_depth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function mod_depth_trials_Callback(hObject, eventdata, handles)
% hObject    handle to mod_depth_trials (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mod_depth_trials as text
%        str2double(get(hObject,'String')) returns contents of mod_depth_trials as a double
[v,handles]=edit_update(handles,hObject,1,50,5,'%d',1);
if v~=handles.trials_depth,
	handles.trials_depth=v;
end;
% Update handles structure
guidata(hObject, handles);
update_trials(handles);

% --- Executes during object creation, after setting all properties.
function mod_depth_trials_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mod_depth_trials (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in previous_set.
function previous_set_Callback(hObject, eventdata, handles)
% hObject    handle to previous_set (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

sel=handles.A_selected;
sel=sel-1;
if sel>=1,
    handles.A_selected=sel;
    set(handles.parameter_set,'String',sprintf('%i',sel));
end;
% Update handles structure
guidata(hObject, handles);
update_plots(handles);


% --- Executes on button press in next_set.
function next_set_Callback(hObject, eventdata, handles)
% hObject    handle to next_set (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
sel=handles.A_selected;
sel=sel+1;
if sel<=handles.successful_trials,
    handles.A_selected=sel;
    set(handles.parameter_set,'String',sprintf('%i',sel));
end;
% Update handles structure
guidata(hObject, handles);
update_plots(handles);


% --- Executes on button press in short_only.
function short_only_Callback(hObject, eventdata, handles)
% hObject    handle to short_only (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

r2=handles.A_r.^2;

sel=1;
best_fom=1e9;
for k=1:handles.successful_trials,
    distr=handles.trial_distr(k,:);
    fom=sum(distr.*r2); % figure of merit
    if fom<best_fom,
        sel=k;
        best_fom=fom;
    end;
end;

handles.A_selected=sel;
set(handles.parameter_set,'String',sprintf('%i',sel));

% Update handles structure
guidata(hObject, handles);
update_plots(handles);



% --- Executes on button press in error_estimate.
function error_estimate_Callback(hObject, eventdata, handles)
% hObject    handle to error_estimate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of error_estimate

% Update handles structure
guidata(hObject, handles);
update_plots(handles);

function update_plots(handles)

rpoi=handles.A_selected;

texp=handles.texp;
vexp=handles.vexp;

bckg=handles.trial_bckg(rpoi,:);
cluster=handles.trial_cluster(rpoi,:);
sim=handles.trial_sim(rpoi,:);
distr=handles.trial_distr(rpoi,:);

% Determine distance constraints for current distribution
rmean=sum(handles.A_r.*distr)/sum(distr);
set(handles.mean_distance,'String',sprintf('%4.1f',rmean));
[m,mpoi]=min(abs(handles.A_r-rmean*ones(size(handles.A_r)))); % mpoi is index of mean distance in r axis
mask=handles.A_r<handles.upper_c_bound*ones(size(handles.A_r));
ndistr=mask.*distr;
ndistr=ndistr/sum(ndistr);
nsum=cumsum(ndistr);
mipoi=1;
stop=0;
k=1;
while ~stop && k<=length(nsum),
    if nsum(k)<(1-handles.constraint_level)/2, 
        mipoi=k;
        k=k+1;
    else
        stop=1;
    end;
end;
mapoi=length(nsum);
stop=0;
k=length(nsum);
while ~stop && k>0,
    if nsum(k)>1-(1-handles.constraint_level)/2, 
        mapoi=k;
        k=k-1;
    else
        stop=1;
    end;
end;

% mipoi=mpoi;
% mapoi=mpoi;
% 
% min_width=length(ndistr);
% for k=1:mpoi-1,
%     for kk=mpoi+1:length(ndistr),
%         if kk-k<min_width,
%             if sum(ndistr(k:kk))>=handles.constraint_level,
%                 mipoi=k;
%                 mapoi=kk;
%                 min_width=kk-k;
%             end;
%         end;
%     end;
% end;

rmin=handles.A_r(mipoi);
rmax=handles.A_r(mapoi);
handles.min_cons=rmin;
handles.max_cons=rmax;
set(handles.constraints,'String',sprintf('%4.1f%s%4.1f',rmin,'-',rmax));

bckg_start=handles.start_vals(rpoi);
set(handles.current_bckg_start,'String',sprintf('%i',bckg_start));
dens=handles.dens_vals(rpoi);
set(handles.current_density,'String',sprintf('%6.3f',dens*handles.calib_density));
depth=handles.depth_vals(rpoi);
set(handles.current_depth,'String',sprintf('%5.3f',depth));
dim=handles.dim_vals(rpoi);
set(handles.current_dim,'String',sprintf('%5.2f',dim));
rmsd=handles.trial_rmsd(rpoi);
set(handles.current_rmsd,'String',sprintf('%8.6f',rmsd));

if handles.copy_flag,
    figure(4); clf;
else
    axes(handles.axes1);
    cla;
end;

plot(texp,real(vexp),'k');
hold on;
plot(handles.A_tdip,bckg,'r');
xlabel('t (탎)');
axis([min(texp),max(texp),min([min(bckg) min(real(vexp))]),max([1.1*max(bckg) 1.1*max(real(vexp))])]);

if handles.copy_flag,
    figure(5); clf;
else
    axes(handles.axes2);
    cla;
end;

plot(handles.A_tdip,cluster,'k');
hold on;
plot(handles.A_tdip,sim,'r');
axis([min(handles.A_tdip),max(handles.A_tdip),min([min(cluster) min(real(sim))]),max([max(cluster) max(real(sim))])]);
xlabel('t (탎)');

if handles.copy_flag,
    figure(6); clf;
else
    axes(handles.axes3);
    cla;
end;

if handles.theor,
    th_distr=handles.th_distr;
    th_r=handles.th_r;
    th_distr=max(distr)*th_distr/max(th_distr);
    plot(handles.th_r,th_distr,'c','LineWidth',1.5);
    hold on;
end;

plot(handles.A_r,distr,'k','LineWidth',1);
hold on
axis([min(handles.A_r),max(handles.A_r),-0.1*max(distr),1.1*max(distr)]);
xlabel('r (nm)');



% plot([rmin rmin],[-0.1*max(distr),1.1*max(distr)],'m:','LineWidth',0.5);
% plot([rmean rmean],[-0.1*max(distr),1.1*max(distr)],'c:','LineWidth',1);
% plot([rmax rmax],[-0.1*max(distr),1.1*max(distr)],'m:','LineWidth',0.5);
% plot([handles.upper_c_bound handles.upper_c_bound],[-0.1*max(distr),1.1*max(distr)],'b:','LineWidth',1);

error_flag=get(handles.error_estimate,'Value'); % show confidence interval for Tikhonov regularization

dlow0=handles.mean_distr-2*handles.distr_std;
dlow=min(handles.trial_distr);
for k=1:length(dlow),
    if dlow(k)<0; dlow(k)=0; end;
end;
dhigh0=handles.mean_distr+2*handles.distr_std;
dhigh=max(handles.trial_distr);

if error_flag && length(dlow)==length(distr) && length(dhigh)==length(distr),
    r=handles.A_r;
    maxihigh=0;
    for k=1:length(r),
        if dhigh(k)>maxihigh, maxihigh=dhigh(k); end;
        if dlow0(k)<0, dlow0(k)=0; end;
        plot([r(k) r(k)],[dlow(k) dhigh(k)],'Color',[0.65 0.65 0.65],'Linewidth',0.5);
    end;
    plot(r,distr,'k','LineWidth',1);
    plot(r,dlow0,'b','LineWidth',1);
    plot(r,dhigh0,'r','LineWidth',1);
    sc=max(dhigh)-min(dlow);
    maxv=max(dhigh)+0.1*sc;
    minv=min(dlow)-0.1*sc;
    axis([min(handles.A_r),max(handles.A_r),minv,maxv]);
end;
if handles.copy_flag,
    set(gca,'FontSize',14);
end;

handles.A_distr=handles.trial_distr(rpoi,:);
handles.A_sim=handles.trial_sim(rpoi,:);
handles.A_cluster=handles.trial_cluster(rpoi,:);
handles.A_dipevo=handles.trial_dipevo(rpoi,:);
handles.A_bckg=handles.trial_bckg(rpoi,:);

handles.copy_flag=0;

% Update handles structure
guidata(gcbo, handles);



% --- Executes on button press in select_best.
function select_best_Callback(hObject, eventdata, handles)
% hObject    handle to select_best (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[mi,sel]=min(handles.trial_rmsd);
handles.A_selected=sel;
set(handles.parameter_set,'String',sprintf('%i',sel));

% Update handles structure
guidata(hObject, handles);
update_plots(handles);




% --- Executes on button press in prune.
function prune_Callback(hObject, eventdata, handles)
% hObject    handle to prune (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~handles.computed, return; end;

tdip=handles.A_tdip;
total_trials=handles.successful_trials;
poi=0;
brmsd=min(handles.trial_rmsd);

trial_bckg=zeros(total_trials,length(tdip));
trial_dipevo=zeros(total_trials,length(tdip));
trial_cluster=zeros(total_trials,length(tdip));
trial_distr=zeros(total_trials,length(handles.A_distr));
trial_sim=zeros(total_trials,length(tdip));
rms_vec=zeros(1,total_trials);
dens_vals=zeros(1,total_trials);
depth_vals=zeros(1,total_trials);
dim_vals=zeros(1,total_trials);
sc_vals=zeros(1,total_trials);

for k=1:total_trials,
    if handles.trial_rmsd(k)<=handles.pruning*brmsd,
        poi=poi+1;
        trial_bckg(poi,:)=handles.trial_bckg(k,:);
        trial_dipevo(poi,:)=handles.trial_dipevo(k,:);
        trial_cluster(poi,:)=handles.trial_cluster(k,:);
        trial_distr(poi,:)=handles.trial_distr(k,:);
        trial_sim(poi,:)=handles.trial_sim(k,:);
        rms_vec(poi)=handles.trial_rmsd(k);
        dens_vals(poi)=handles.dens_vals(k);
        depth_vals(poi)=handles.depth_vals(k);
        dim_vals(poi)=handles.dim_vals(k);
        sc_vals(poi)=handles.sc_vals(k);
    end;
end;

handles.successful_trials=poi;
trial_bckg=trial_bckg(1:poi,:);
trial_dipevo=trial_dipevo(1:poi,:);
trial_cluster=trial_cluster(1:poi,:);
trial_distr=trial_distr(1:poi,:);
trial_sim=trial_sim(1:poi,:);
rms_vec=rms_vec(1:poi);
handles.dens_vals=dens_vals(1:poi);
handles.depth_vals=depth_vals(1:poi);
handles.dim_vals=dim_vals(1:poi);
handles.sc_vals=sc_vals(1:poi);

[min_rms,rpoi]=min(rms_vec);
handles.dens=dens_vals(rpoi);
handles.mod_depth=depth_vals(rpoi);
distr_std=std(trial_distr);
mean_distr=mean(trial_distr);
handles.mean_distr=mean_distr;
handles.distr_std=distr_std;
handles.A_selected=rpoi;

handles.trial_distr=trial_distr;
handles.trial_sim=trial_sim;
handles.trial_cluster=trial_cluster;
handles.trial_dipevo=trial_dipevo;
handles.trial_bckg=trial_bckg;
handles.trial_rmsd=rms_vec;

set(handles.param_line,'String',sprintf('%i%s%i%s%5.2f%s',poi,' out of ',total_trials,' data sets within ',handles.pruning,' times best r.m.s.d.'));
set(handles.parameter_set,'String',sprintf('%i',rpoi));

% Update handles structure
guidata(hObject, handles);
update_plots(handles);


% --- Executes on button press in rough_grid.
function rough_grid_Callback(hObject, eventdata, handles)
% hObject    handle to rough_grid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global main_handles

% load validation_data
defaults=load('validation_defaults.dat');

if handles.computed,
    rpoi=handles.A_selected;
    dens0=handles.dens_vals(rpoi);
    depth0=handles.depth_vals(rpoi);
    dim0=handles.dim_vals(rpoi);
else,
    dens0=handles.dens;
    depth0=handles.mod_depth;
    dim0=handles.hom_dim;
end;

if get(handles.noise_test,'Value'),
    handles.trials_noise=5;
    handles.factor_noise=2.0;
    set(handles.noise_factor,'String',sprintf('%4.1f',handles.factor_noise));
else,
    handles.trials_noise=1;
    handles.factor_noise=1.0;
    set(handles.noise_factor,'String',sprintf('%4.1f',handles.factor_noise));
end;
set(handles.noise_trials,'String',num2str(handles.trials_noise));
if get(handles.background_density,'Value'),
    handles.trials_bckg=11;
    handles.min_bckg=0.5*dens0;
    pstr=sprintf('%6.4f',handles.min_bckg*handles.calib_density);
    set(handles.bckg_min,'String',pstr);
    handles.max_bckg=2*dens0;
    pstr=sprintf('%6.4f',handles.max_bckg*handles.calib_density);
    set(handles.bckg_max,'String',pstr);
else,
    handles.trials_bckg=1;
    handles.min_bckg=dens0;
    pstr=sprintf('%6.4f',handles.min_bckg*handles.calib_density);
    set(handles.bckg_min,'String',pstr);
    handles.max_bckg=dens0;
    pstr=sprintf('%6.4f',handles.max_bckg);
    set(handles.bckg_max,'String',pstr);
end;
set(handles.bckg_trials,'String',num2str(handles.trials_bckg));
if get(handles.vary_depth,'Value'),
    handles.trials_depth=11;
    handles.mod_depth_min=0.5*depth0;
    pstr=sprintf('%5.3f',handles.mod_depth_min);
    set(handles.min_mod_depth,'String',pstr);
    handles.mod_depth_max=2*depth0;
    pstr=sprintf('%5.3f',handles.mod_depth_max);
    set(handles.max_mod_depth,'String',pstr);
else,
    handles.trials_depth=1;
    handles.mod_depth_min=depth0;
    pstr=sprintf('%5.3f',handles.mod_depth_min);
    set(handles.min_mod_depth,'String',pstr);
    handles.mod_depth_max=depth0;
    pstr=sprintf('%5.3f',handles.mod_depth_max);
    set(handles.max_mod_depth,'String',pstr);
end;
set(handles.mod_depth_trials,'String',num2str(handles.trials_depth));
if get(handles.background_dim,'Value'),
    handles.trials_dim=6;
    handles.min_dim=2;
    pstr=sprintf('%5.2f',handles.min_dim);
    set(handles.dim_min,'String',pstr);
    handles.max_dim=3;
    pstr=sprintf('%5.2f',handles.max_dim);
    set(handles.dim_max,'String',pstr);
else
    handles.trials_dim=1;
    handles.min_dim=dim0;
    pstr=sprintf('%5.2f',handles.min_dim);
    set(handles.dim_min,'String',pstr);
    handles.max_dim=dim0;
    pstr=sprintf('%5.2f',handles.max_dim);
    set(handles.dim_max,'String',pstr);
end;
set(handles.dim_trials,'String',num2str(handles.trials_dim));
if get(handles.checkbox_bckg_start,'Value'),
    handles.trials_bckg_start=6;
    handles.min_bckg_start=defaults(18);
    pstr=sprintf('%i',handles.min_bckg_start);
    set(handles.edit_bckg_start_min,'String',pstr);
    handles.max_bckg_start=defaults(19);
    pstr=sprintf('%i',handles.max_bckg_start);
    set(handles.edit_bckg_start_max,'String',pstr);
else
    handles.trials_bckg_start=1;
    handles.min_bckg_start=main_handles.bckg_start;
    pstr=sprintf('%i',handles.min_bckg_start);
    set(handles.edit_bckg_start_min,'String',pstr);
    handles.max_bckg_start=main_handles.bckg_start;
    pstr=sprintf('%i',handles.max_bckg_start);
    set(handles.edit_bckg_start_max,'String',pstr);
end;
set(handles.edit_bckg_start_trials,'String',num2str(handles.trials_bckg_start));
handles.trials_total=handles.trials_noise*handles.trials_bckg*handles.trials_depth*handles.trials_dim*handles.trials_bckg_start;
set(handles.total_trials,'String',num2str(handles.trials_total));
handles.computed=0;
set(handles.prune,'Enable','off');
set(handles.status_line,'String',sprintf('%s%6.4f%s%5.3f%s%5.2f','Initial density: ',handles.dens*handles.calib_density,', Initial depth: ',handles.mod_depth,', Initial dimension: ',handles.hom_dim));

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in fine_grid.
function fine_grid_Callback(hObject, eventdata, handles)
% hObject    handle to fine_grid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global main_handles
% load validation_data
defaults=load('validation_defaults.dat');

if handles.computed,
    rpoi=handles.A_selected;
    dens0=handles.dens_vals(rpoi);
    dens_min=min(handles.dens_vals);
    dens_max=max(handles.dens_vals);
    depth0=handles.depth_vals(rpoi);
    depth_min=min(handles.depth_vals);
    depth_max=max(handles.depth_vals);
    dim0=handles.dim_vals(rpoi);
    dim_min=min(handles.dim_vals);
    dim_max=max(handles.dim_vals);
else,
    dens0=handles.dens;
    dens_min=0.8*handles.dens;
    dens_max=1.2*handles.dens;
    depth0=handles.mod_depth;
    depth_min=0.8*handles.mod_depth;
    depth_max=1.2*handles.mod_depth;
    dim0=handles.hom_dim;
    dim_min=dim0-0.2;
    dim_max=dim0+0.2;
end;

if get(handles.noise_test,'Value'),
    handles.trials_noise=10;
    handles.factor_noise=2.0;
    set(handles.noise_factor,'String',sprintf('%4.1f',handles.factor_noise));
else,
    handles.trials_noise=1;
    handles.factor_noise=1.0;
    set(handles.noise_factor,'String',sprintf('%4.1f',handles.factor_noise));
end;
set(handles.noise_trials,'String',num2str(handles.trials_noise));
if get(handles.background_density,'Value'),
    handles.trials_bckg=11;
    handles.min_bckg=dens_min;
    pstr=sprintf('%6.4f',handles.min_bckg*handles.calib_density);
    set(handles.bckg_min,'String',pstr);
    handles.max_bckg=dens_max;
    pstr=sprintf('%6.4f',handles.max_bckg*handles.calib_density);
    set(handles.bckg_max,'String',pstr);
else,
    handles.trials_bckg=1;
    handles.min_bckg=dens0;
    pstr=sprintf('%6.4f',handles.min_bckg*handles.calib_density);
    set(handles.bckg_min,'String',pstr);
    handles.max_bckg=dens0;
    pstr=sprintf('%6.4f',handles.max_bckg*handles.calib_density);
    set(handles.bckg_max,'String',pstr);
end;
set(handles.bckg_trials,'String',num2str(handles.trials_bckg));
if get(handles.vary_depth,'Value'),
    handles.trials_depth=11;
    handles.mod_depth_min=depth_min;
    pstr=sprintf('%5.3f',handles.mod_depth_min);
    set(handles.min_mod_depth,'String',pstr);
    handles.mod_depth_max=depth_max;
    pstr=sprintf('%5.3f',handles.mod_depth_max);
    set(handles.max_mod_depth,'String',pstr);
else,
    handles.trials_depth=1;
    handles.mod_depth_min=depth0;
    pstr=sprintf('%5.3f',handles.mod_depth_min);
    set(handles.min_mod_depth,'String',pstr);
    handles.mod_depth_max=depth0;
    pstr=sprintf('%5.3f',handles.mod_depth_max);
    set(handles.max_mod_depth,'String',pstr);
end;
set(handles.mod_depth_trials,'String',num2str(handles.trials_depth));
if get(handles.background_dim,'Value'),
    handles.trials_dim=9;
    handles.min_dim=dim0-0.2;
    pstr=sprintf('%5.2f',handles.min_dim);
    set(handles.dim_min,'String',pstr);
    handles.max_dim=dim0+0.2;
    pstr=sprintf('%5.2f',handles.max_dim);
    set(handles.dim_max,'String',pstr);
else,
    handles.trials_dim=1;
    handles.min_dim=handles.hom_dim;
    pstr=sprintf('%5.2f',handles.min_dim);
    set(handles.dim_min,'String',pstr);
    handles.max_dim=handles.hom_dim;
    pstr=sprintf('%5.2f',handles.max_dim);
    set(handles.dim_max,'String',pstr);
end;
set(handles.dim_trials,'String',num2str(handles.trials_dim));
if get(handles.checkbox_bckg_start,'Value'),
    handles.trials_bckg_start=11;
    handles.min_bckg_start=defaults(18);
    pstr=sprintf('%i',handles.min_bckg_start);
    set(handles.edit_bckg_start_min,'String',pstr);
    handles.max_bckg_start=defaults(19);
    pstr=sprintf('%i',handles.max_bckg_start);
    set(handles.edit_bckg_start_max,'String',pstr);
else
    handles.trials_bckg_start=1;
    handles.min_bckg_start=main_handles.bckg_start;
    pstr=sprintf('%i',handles.min_bckg_start);
    set(handles.edit_bckg_start_min,'String',pstr);
    handles.max_bckg_start=main_handles.bckg_start;
    pstr=sprintf('%i',handles.max_bckg_start);
    set(handles.edit_bckg_start_max,'String',pstr);
end;
set(handles.edit_bckg_start_trials,'String',num2str(handles.trials_bckg_start));
handles.trials_total=handles.trials_noise*handles.trials_bckg*handles.trials_depth*handles.trials_dim*handles.trials_bckg_start;
set(handles.total_trials,'String',num2str(handles.trials_total));
handles.computed=0;
set(handles.prune,'Enable','off');
set(handles.status_line,'String',sprintf('%s%6.4f%s%5.3f%s%5.2f','Initial density: ',handles.dens*handles.calib_density,', Initial depth: ',handles.mod_depth,', Initial dimension: ',handles.hom_dim));

% Update handles structure
guidata(hObject, handles);



function prune_level_Callback(hObject, eventdata, handles)
% hObject    handle to prune_level (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of prune_level as text
%        str2double(get(hObject,'String')) returns contents of prune_level as a double
[v,handles]=edit_update(handles,hObject,1.01,3,1.5,'%5.2f',0);
if v~=handles.pruning,
	handles.pruning=v;
end;
% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function prune_level_CreateFcn(hObject, eventdata, handles)
% hObject    handle to prune_level (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in statusbar_off.
function statusbar_off_Callback(hObject, eventdata, handles)
% hObject    handle to statusbar_off (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of statusbar_off




% --- Executes on button press in compact.
function compact_Callback(hObject, eventdata, handles)
% hObject    handle to compact (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

r=handles.A_r;

sel=1;
best_mom2=1e9;
for k=1:handles.successful_trials,
    distr=handles.trial_distr(k,:);
    mom1=sum(distr.*r)/sum(distr);
    r2=r-mom1*ones(size(r));
    mom2=sum(distr.*r2.^2)/sum(distr); % figure of merit
    if mom2<best_mom2,
        sel=k;
        best_mom2=mom2;
    end;
end;

handles.A_selected=sel;
set(handles.parameter_set,'String',sprintf('%i',sel));

% Update handles structure
guidata(hObject, handles);
update_plots(handles);




% --- Executes on button press in copy.
function copy_Callback(hObject, eventdata, handles)
% hObject    handle to copy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.copy_flag=1;
% Update handles structure
guidata(hObject, handles);
update_plots(handles);

function parameter_set_Callback(hObject, eventdata, handles)
% hObject    handle to bckg_trials (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of bckg_trials as text
%        str2double(get(hObject,'String')) returns contents of bckg_trials as a double
[v,handles]=edit_update(handles,hObject,1,handles.successful_trials,handles.A_selected,'%d',1);
if v~=handles.A_selected,
	handles.A_selected=v;
end;
% Update handles structure
guidata(hObject, handles);
update_plots(handles);


% --- Executes during object creation, after setting all properties.
function parameter_set_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bckg_trials (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on button press in worst.
function worst_Callback(hObject, eventdata, handles)
% hObject    handle to worst (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[ma,sel]=max(handles.trial_rmsd);
handles.A_selected=sel;
set(handles.parameter_set,'String',sprintf('%i',sel));

% Update handles structure
guidata(hObject, handles);
update_plots(handles);




function upper_bound_Callback(hObject, eventdata, handles)
% hObject    handle to upper_bound (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of upper_bound as text
%        str2double(get(hObject,'String')) returns contents of upper_bound as a double
[v,handles]=edit_update(handles,hObject,min(handles.A_r),max(handles.A_r),handles.upper_c_bound,'%5.2f',0);
if v~=handles.upper_c_bound,
	handles.upper_c_bound=v;
end;
% Update handles structure
guidata(hObject, handles);
update_plots(handles);


% --- Executes during object creation, after setting all properties.
function upper_bound_CreateFcn(hObject, eventdata, handles)
% hObject    handle to upper_bound (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in checkbox_bckg_start.
function checkbox_bckg_start_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_bckg_start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_bckg_start

guidata(hObject, handles);
update_trials(handles);


function edit_bckg_start_trials_Callback(hObject, eventdata, handles)
% hObject    handle to edit_bckg_start_trials (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_bckg_start_trials as text
%        str2double(get(hObject,'String')) returns contents of edit_bckg_start_trials as a double

[v,handles]=edit_update(handles,hObject,1,50,5,'%d',1);
handles.trials_bckg_start=v;
% Update handles structure
guidata(hObject, handles);
update_trials(handles);

% --- Executes during object creation, after setting all properties.
function edit_bckg_start_trials_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_bckg_start_trials (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_bckg_start_min_Callback(hObject, eventdata, handles)
% hObject    handle to edit_bckg_start_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_bckg_start_min as text
%        str2double(get(hObject,'String')) returns contents of edit_bckg_start_min as a double

global main_handles

% load validation_data
[v,handles]=edit_update(handles,hObject,0,main_handles.cutoff,0,'%d',1);
handles.min_bckg_start=v;
% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit_bckg_start_min_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_bckg_start_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_bckg_start_max_Callback(hObject, eventdata, handles)
% hObject    handle to edit_bckg_start_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_bckg_start_max as text
%        str2double(get(hObject,'String')) returns contents of edit_bckg_start_max as a double

global main_handles
% load validation_data
[v,handles]=edit_update(handles,hObject,0,main_handles.cutoff,0,'%d',1);
handles.max_bckg_start=v;
% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function edit_bckg_start_max_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_bckg_start_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
