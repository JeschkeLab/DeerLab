function handback=set_defaults(handles)
% Sets all internal status variable to their default values

% Control panel: Original data
set(handles.dual_display,'Value',0);
set(handles.imaginary,'Value',1);
set(handles.imaginary,'String','imaginary');

handles.zerotime=0;
pstr=sprintf('%d',handles.zerotime);
set(handles.zt_edit,'String',pstr);

handles.phase=0;
pstr=sprintf('%d',handles.phase);
set(handles.phase_edit,'String',pstr);

handles.bckg_start=round(max(handles.A_texp)/4);
pstr=sprintf('%d',handles.bckg_start);
set(handles.bckg_edit,'String',pstr);

handles.cutoff=max(handles.A_texp);
handles.cutoff_suggestion=max(handles.A_texp);
pstr=sprintf('%d',handles.cutoff);
set(handles.cutoff_edit,'String',pstr);

% Control panel: Dipolar evolution
set(handles.dip_time_domain,'Value',1);
set(handles.dip_spectrum,'Value',0);
set(handles.long_pass,'Value',0);

handles.longpass_min=1.6;
pstr=sprintf('%4.1f',handles.longpass_min);
set(handles.long_pass_edit,'String',pstr);

handles.zoom=1.0;
pstr=sprintf('%6.3f',handles.zoom);
set(handles.zoom_edit,'String',pstr);

set(handles.exci_bandwidth_corr,'Value',0);
handles.bandwidth=16;
pstr=sprintf('%6.1f',handles.bandwidth);
set(handles.exci_bandwidth_edit,'String',pstr);

set(handles.checkbox_ghost,'Value',0);
handles.spins_per_object=3;
set(handles.edit_oligomerization,'String',sprintf('%d',3));

% Control panel: Distance distribution
handles.rmin=1.5;
pstr=sprintf('%4.1f',handles.rmin);
set(handles.distr_start_edit,'String',pstr);

handles.rmax=8.0;
pstr=sprintf('%4.1f',handles.rmax);
set(handles.distr_end_edit,'String',pstr);

set(handles.L_curve,'Value',0);

handles.regpar=1.0;
set(handles.regpar_edit,'String',num2str(handles.regpar,handles.regpar_edit_strformat));

handles.regpar_sel=3;

% Control panel: Data sets, not automatic defaults!

% Control panel: Background model
handles.bckg_rms_value=1.0e6;
set(handles.bckg_rms,'String','n.a.');

handles.hom_density=1.0;
set(handles.bckg_density,'String','n.a.');

set(handles.bckg_none,'Value',0);
set(handles.bckg_homogeneous,'Value',1);
handles.hom_dim=3;
pstr=sprintf('%d',handles.hom_dim);
set(handles.bckg_dim_edit,'String',pstr);
set(handles.bckg_fit_dim,'Value',0);
set(handles.bckg_poly,'Value',0);
set(handles.bckg_poly,'Value',0);
handles.poly_order=5;
pstr=sprintf('%d',handles.poly_order);
set(handles.bckg_poly_order,'String',pstr);
set(handles.bckg_exp,'Value',0);
set(handles.radiobutton_deernet_bckg,'Value',0);
set(handles.checkbox_deernet_error_bckg,'Value',0);
set(handles.radiobutton_deernet_bckg,'Enable','off');
set(handles.checkbox_deernet_error_bckg,'Enable','off');

% Control panel: Distance analysis
handles.fit_rms_value=1.0e6;
set(handles.distr_rms,'String','n.a.');
handles.r_MA=5;
set(handles.r_mean,'String','n.a.');
handles.sr_MA=0.5;
set(handles.sr,'String','n.a.');
handles.n_spins=2.0;
set(handles.num_spins,'String','n.a.');
handles.moddepth_suppression=1;

set(handles.select_APT,'Value',1);
handles.DDS=0.2;
pstr=sprintf('%0.2g',handles.DDS);
set(handles.DDS_filter,'String',pstr);

set(handles.select_deernet,'Value',0);
set(handles.select_Tikhonov,'Value',0);
set(handles.select_L_curve,'Value',1);
handles.regpars=1;

set(handles.select_model,'Value',0);
set(handles.user_model_list,'Value',1);

handles.user_fit_rmin=1.5;
handles.user_fit_rmax=10;
handles.user_fit_dr=0.02;

handles.A_low = [];
handles.A_high= [];

% initialize DEERNet
handles.A_deernet_t = [];
handles.A_deernet_sim = [];
handles.A_deernet_ff = [];
handles.A_deernet_bckg = [];
handles.A_deernet_vexp = [];

handles.deernet_ensemble_sim = [];
handles.deernet_ensemble_ff = [];
handles.deernet_ensemble_bckg = [];


% prepare for comparative mode
handles.A_curr_r = [];
handles.A_curr_distr = [];
handles.A_prev_r = [];
handles.A_prev_distr = [];
handles.new_distr = 1;

handles.best_rmsd = 1e6;
handles.user_model_trial = false;

handback=read_user_model(handles);

handback.A_prev_bckg_mode = '';
handback.A_prev_bckg_details = '';
handback.A_prev_bckg = [];
handback.A_prev_bckg_t = [];

handback.A_curr_bckg_mode = '';
handback.A_curr_bckg_details = '';
handback.A_curr_bckg = [];
handback.A_curr_bckg_t = [];




