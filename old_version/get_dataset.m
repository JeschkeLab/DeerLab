function handback=get_dataset(handles)
%
% Load data set in the selected format and save it in data structure handles
%

handback=handles;

% Store current data set as data set B
handles.theor=0;
handles.B_texp=handles.A_texp;
handles.B_vexp=handles.A_vexp;
handles.B_vb=handles.A_vb;
handles.B_bas_name=handles.bas_name;
handles.B_tplot=handles.texp;
handles.B_vplot=handles.vexp;
handles.B_tdip=handles.A_tdip;
handles.B_dipevo=handles.A_dipevo;
handles.B_cluster=handles.A_cluster;
handles.B_spc=handles.A_spc;
handles.B_ny=handles.A_ny;
handles.B_r=handles.A_r;
handles.B_distr=handles.A_distr;
handles.B_depth=handles.A_depth;
handles.updated=0;
handles.validation_mode=0;
set(handles.validate_Tikhonov,'Enable','off');
set(handles.manual_bckg,'Value',0);
handles.locked_loaded=false;

flag=get(handles.reset_on_load,'Value');
if flag, handles=set_defaults(handles); end % set all status variables to defaults

data_format=get(handles.format_winepr,'Value'); % 0 for Elexsys or ASCII, 1 for WINEPR
flag=get(handles.format_ascii,'Value');
if flag, data_format=2; end
if get(handles.radiobutton_DeerAnalysis,'Value')
    data_format=3;
end
switch data_format
    case 0 % Elexsys
        ext='*.DTA;';
    case 1 % WINEPR
        ext='*.spc;';
    case 2 % ASCII
        ext='*.txt;*.asc;*.dat';
    case 3 % DeerAnalysis
        ext='*.txt;*.dat';
end
[fname,pname]=uigetfile(ext,'Load experimental dataset');
if isequal(fname,0)||isequal(pname,0)
    handback=handles;
    return; 
end
cd(pname);

% separate filename from extension
dots=findstr('.',fname);
ext_pos=length(fname);
if ~isempty(dots), ext_pos=dots(length(dots))-1; end
bas_name=fname(1:ext_pos); % name without extension
extension=fname(ext_pos+2:length(fname)); % only extension
pfname=[pname fname]; % pathname
pbname=[pname bas_name]; % pathname without extension
thname=[bas_name '_theor.distr']; % name of file with theoretical distribution
header_lines=0;
handles.project_dir=pname;
handles.bas_name=bas_name;
handles.source_file=[pname fname];

success=true;
switch data_format
    case 0 % Elexsys
        [x,y,z,vb]=get_elexsys(bas_name);
    case 1 % WINEPR
        [x,y,z,vb]=get_WINEPR(bas_name);
    case 2 % ASCII
        dset=load(fname,'-ascii');
        [mtest,ntest]=size(dset); % determine size of loaded data array
        x=dset(:,handles.ASCII_t_column);
        x=x';
        z=dset(:,handles.ASCII_real_column);
        if ntest>=handles.ASCII_imag_column
            z=z+1i*dset(:,handles.ASCII_imag_column);
        end
        z=real(z')-1i*imag(z');
        %z=z';
        y=[];
        vb=[];
        % check, if test data set and, if so, load and store theoretical
        % distribution
        if exist(thname,'file')
            % disp('Test data set');
            dset2=load(thname,'-ascii');
            handles.th_r=dset2(:,1);
            handles.th_distr=dset2(:,2);
            handles.theor=1;
        else
            handles.theor=0;
        end
    case 3
        % determine true basis name
        success=false;
        pos=strfind(bas_name,'_');
        if ~isempty(pos)
            bas_name=bas_name(1:pos(end)-1);
        end
        % check if all files exist
        success=true;
        bckg_file=strcat(bas_name,'_bckg.dat');
        if ~exist(bckg_file,'file')
            success=false;
        end
        fit_file=strcat(bas_name,'_fit.dat');
        if ~exist(fit_file,'file')
            success=false;
        end
        distr_file=strcat(bas_name,'_distr.dat');
        if ~exist(distr_file,'file')
            success=false;
        end
        if success
            dset=load(bckg_file,'-ascii');
            x=dset(:,1);
            x=x';
            z=dset(:,2);
            z=z';
            handles.loaded_bckg=dset(:,3)';
            dset=load(fit_file,'-ascii');
            handles.loaded_tdip=dset(:,1)';
            handles.loaded_ff=dset(:,2)';
            handles.loaded_sim=dset(:,3)';
            dset=load(distr_file,'-ascii');
            [~,nd]=size(dset);
            handles.loaded_rax=dset(:,1)';
            handles.loaded_distr=dset(:,2)';
            if nd>3
                handles.loaded_distr_lower=dset(:,3)';
                handles.loaded_distr_upper=dset(:,4)';
            end
            vb=[];
        end
end
if ~success
    set(handles.status_line,'String','Complete set of DeerAnalysis output files not available.');
    return
end
figname=['DeerAnalysis 2018 - ' fname]; % tell user, which file is current
set(handles.main_figure,'Name',figname);
dx=x(2)-x(1);
% correct usec time axis to ns, if required
if dx<0.1
    x=1000*x;
end
if ~get(handles.checkbox_no_analysis,'Value')
    x=x-x(1)*ones(size(x));
else
    handles.locked_loaded=true;
end
[mtest,~]=size(z); % determine size of loaded data array
if mtest==2
    handles.ctvt=1;
    dtype='Variable-time DEER. ';
else
    handles.ctvt=0;
    dtype='Constant-time DEER. ';
end
cflag=round(1000*sum(abs(imag(z)))/sum(abs(real(z)))); % Determine, if imaginary part is significant
if cflag
    handles.cmplx=1;
    dform=' complex ';
else
    handles.cmplx=0;
    dform=' real ';
end
nexp=length(x);
% msg=sprintf('%s%d%s%s',dtype,nexp,dform,'data points.');
handles.t_orig=x;
handles.v_orig=z;

[texp,vexp,zt,phi,imo,dt]=pre_process(x,z);
handles.dt=dt;
handles.min_dt=dt;
nmax=length(texp);
handles.max_dt=floor(nmax/64)*dt;

set(handles.data_set_B,'String',handles.B_bas_name);
handles.bas_name=bas_name;
set(handles.data_set_A,'String',bas_name);
currdir=pwd;

handles.A_texp=texp;
handles.A_vexp=vexp;
handles.A_vb=vb;

reset_flag=get(handles.reset_on_load,'Value');
if reset_flag
    handles=set_defaults(handles);
end

% Update handles structure
guidata(handles.main_figure, handles);

zt_string=sprintf('%d',zt);
set(handles.zt_edit,'String',zt_string);
% zt=get_zerotime(handles,x,real(z));
handles.zerotime=zt;

phaseflag=get(handles.autophase,'Value');
if ~phaseflag
    vexp=vexp/max(real(vexp));
    handles.A_vexp=vexp;
    phi=0;
    imo=0;
end
handles.phase=phi;
handles.imo=imo;
pstr=sprintf('%6.1f',phi*180/pi);
set(handles.phase_edit,'String',pstr);
    
handles=update_kernel(handles,nmax);

% Empty result variables
handles.APT=[];
handles.r_APT=[];

% Update handles structure
guidata(handles.main_figure, handles);

handback=handles;

