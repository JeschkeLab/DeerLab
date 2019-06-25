function handback=read_user_model(handles),
%
% Reads the parameter list of the currently selected user model
%

% Disable all parameter controls
set(handles.sel_par1, 'Value',0);
set(handles.sel_par1, 'String','Par. 1');
set(handles.sel_par1, 'Enable','off');
set(handles.par1_edit,'ForegroundColor','r');
set(handles.par1_edit,'String','n.a.');
set(handles.par1_edit,'Enable','off');
set(handles.sel_par2, 'Value',0);
set(handles.sel_par2, 'String','Par. 2');
set(handles.sel_par2, 'Enable','off');
set(handles.par2_edit,'ForegroundColor','r');
set(handles.par2_edit,'String','n.a.');
set(handles.par2_edit,'Enable','off');
set(handles.sel_par3, 'Value',0);
set(handles.sel_par3, 'String','Par. 3');
set(handles.sel_par3, 'Enable','off');
set(handles.par3_edit,'ForegroundColor','r');
set(handles.par3_edit,'String','n.a.');
set(handles.par3_edit,'Enable','off');
set(handles.sel_par4, 'Value',0);
set(handles.sel_par4, 'String','Par. 4');
set(handles.sel_par4, 'Enable','off');
set(handles.par4_edit,'ForegroundColor','r');
set(handles.par4_edit,'String','n.a.');
set(handles.par4_edit,'Enable','off');
set(handles.sel_par5, 'Value',0);
set(handles.sel_par5, 'String','Par. 5');
set(handles.sel_par5, 'Enable','off');
set(handles.par5_edit,'ForegroundColor','r');
set(handles.par5_edit,'String','n.a.');
set(handles.par5_edit,'Enable','off');
set(handles.sel_par6, 'Value',0);
set(handles.sel_par6, 'String','Par. 6');
set(handles.sel_par6, 'Enable','off');
set(handles.par6_edit,'ForegroundColor','r');
set(handles.par6_edit,'String','n.a.');
set(handles.par6_edit,'Enable','off');
set(handles.sel_par7, 'Value',0);
set(handles.sel_par7, 'String','Par. 7');
set(handles.sel_par7, 'Enable','off');
set(handles.par7_edit,'ForegroundColor','r');
set(handles.par7_edit,'String','n.a.');
set(handles.par7_edit,'Enable','off');
set(handles.sel_par8, 'Value',0);
set(handles.sel_par8, 'String','Par. 8');
set(handles.sel_par8, 'Enable','off');
set(handles.par8_edit,'ForegroundColor','r');
set(handles.par8_edit,'String','n.a.');
set(handles.par8_edit,'Enable','off');

handles.model_pars_labels=cell(1,8);
handles.model_pars=zeros(1,8);
handles.model_mask=zeros(1,8);
handles.model_defaults=zeros(1,8);
handles.model_lower_bounds=zeros(1,8);
handles.model_upper_bounds=zeros(1,8);
handles.model_mode=0;
enabled=8; % number of parameters that are enabled by default

model_list=get(handles.user_model_list,'String');
selection=get(handles.user_model_list,'Value');
model=deblank(model_list(selection,:));
fname=[handles.model_path model '.m'];
handles.user_model=model;
fid=fopen(fname,'rt');
while feof(fid) ==0,
    tline=fgetl(fid);
    % Check for number of enabled fit parameters
    t1=findstr(tline,'#enable#');
    if ~isempty(t1),
        if tline(1)=='%',
            defline=tline(2:end);
            [mnemo,defline]=strtok(defline);
            [number,defline]=strtok(defline);
            enabled=str2double(number);
        end;
    end;
    % Check for parameter 1
    t1=findstr(tline,'par(1)');
    if ~isempty(t1),
        if tline(1)=='%',
            defline=tline(t1(1)+6:length(tline));
            [mnemo,defline]=strtok(defline);
            handles.model_pars_labels{1}=mnemo;
            [default,defline]=strtok(defline);
            handles.model_defaults(1)=str2double(default);
            [lobound,defline]=strtok(defline);
            handles.model_lower_bounds(1)=str2double(lobound);
            [hibound,defline]=strtok(defline);
            handles.model_upper_bounds(1)=str2double(hibound);
			set(handles.sel_par1,'String',mnemo);
            if enabled>=1,
            	set(handles.sel_par1,'Value',1);
            else
            	set(handles.sel_par1,'Value',0);
            end;
            set(handles.sel_par1,'Enable','on');
			set(handles.sel_par1,'TooltipString',defline);
			set(handles.par1_edit,'String',default);
			set(handles.par1_edit,'Enable','on');
			set(handles.par1_edit,'TooltipString',defline);
        end;
    end;
    % Check for parameter 2
    t1=findstr(tline,'par(2)');
    if ~isempty(t1),
        if tline(1)=='%',
            defline=tline(t1(1)+6:length(tline));
            [mnemo,defline]=strtok(defline);
            handles.model_pars_labels{2}=mnemo;
            [default,defline]=strtok(defline);
            handles.model_defaults(2)=str2double(default);
            [lobound,defline]=strtok(defline);
            handles.model_lower_bounds(2)=str2double(lobound);
            [hibound,defline]=strtok(defline);
            handles.model_upper_bounds(2)=str2double(hibound);
			set(handles.sel_par2,'String',mnemo);
            if enabled>=2,
            	set(handles.sel_par2,'Value',1);
            else
            	set(handles.sel_par2,'Value',0);
            end;
			set(handles.sel_par2,'TooltipString',defline);
			set(handles.sel_par2,'Enable','on');
			set(handles.par2_edit,'String',default);
			set(handles.par2_edit,'Enable','on');
			set(handles.par2_edit,'TooltipString',defline);
        end;
    end;
    % Check for parameter 3
    t1=findstr(tline,'par(3)');
    if ~isempty(t1),
        if tline(1)=='%',
            defline=tline(t1(1)+6:length(tline));
            [mnemo,defline]=strtok(defline);
            handles.model_pars_labels{3}=mnemo;
            [default,defline]=strtok(defline);
            handles.model_defaults(3)=str2double(default);
            [lobound,defline]=strtok(defline);
            handles.model_lower_bounds(3)=str2double(lobound);
            [hibound,defline]=strtok(defline);
            handles.model_upper_bounds(3)=str2double(hibound);
			set(handles.sel_par3,'String',mnemo);
            if enabled>=3,
            	set(handles.sel_par3,'Value',1);
            else
            	set(handles.sel_par3,'Value',0);
            end;
			set(handles.sel_par3,'Enable','on');
			set(handles.sel_par3,'TooltipString',defline);
			set(handles.par3_edit,'String',default);
			set(handles.par3_edit,'Enable','on');
			set(handles.par3_edit,'TooltipString',defline);
        end;
    end;
    % Check for parameter 4
    t1=findstr(tline,'par(4)');
    if ~isempty(t1),
        if tline(1)=='%',
            defline=tline(t1(1)+6:length(tline));
            [mnemo,defline]=strtok(defline);
            handles.model_pars_labels{4}=mnemo;
            [default,defline]=strtok(defline);
            handles.model_defaults(4)=str2double(default);
            [lobound,defline]=strtok(defline);
            handles.model_lower_bounds(4)=str2double(lobound);
            [hibound,defline]=strtok(defline);
            handles.model_upper_bounds(4)=str2double(hibound);
			set(handles.sel_par4,'String',mnemo);
            if enabled>4,
            	set(handles.sel_par4,'Value',1);
            else
            	set(handles.sel_par4,'Value',0);
            end;
			set(handles.sel_par4,'Enable','on');
			set(handles.sel_par4,'TooltipString',defline);
			set(handles.par4_edit,'String',default);
			set(handles.par4_edit,'Enable','on');
			set(handles.par4_edit,'TooltipString',defline);
        end;
    end;
    % Check for parameter 5
    t1=findstr(tline,'par(5)');
    if ~isempty(t1),
        if tline(1)=='%',
            defline=tline(t1(1)+6:length(tline));
            [mnemo,defline]=strtok(defline);
            handles.model_pars_labels{5}=mnemo;
            [default,defline]=strtok(defline);
            handles.model_defaults(5)=str2double(default);
            [lobound,defline]=strtok(defline);
            handles.model_lower_bounds(5)=str2double(lobound);
            [hibound,defline]=strtok(defline);
            handles.model_upper_bounds(5)=str2double(hibound);
			set(handles.sel_par5,'String',mnemo);
            if enabled>=5,
            	set(handles.sel_par5,'Value',1);
            else
            	set(handles.sel_par5,'Value',0);
            end;
			set(handles.sel_par5,'Enable','on');
			set(handles.sel_par5,'TooltipString',defline);
			set(handles.par5_edit,'String',default);
			set(handles.par5_edit,'Enable','on');
			set(handles.par5_edit,'TooltipString',defline);
        end;
    end;
    % Check for parameter 6
    t1=findstr(tline,'par(6)');
    if ~isempty(t1),
        if tline(1)=='%',
            defline=tline(t1(1)+6:length(tline));
            [mnemo,defline]=strtok(defline);
            handles.model_pars_labels{6}=mnemo;
            [default,defline]=strtok(defline);
            handles.model_defaults(6)=str2double(default);
            [lobound,defline]=strtok(defline);
            handles.model_lower_bounds(6)=str2double(lobound);
            [hibound,defline]=strtok(defline);
            handles.model_upper_bounds(6)=str2double(hibound);
			set(handles.sel_par6,'String',mnemo);
            if enabled>=6,
            	set(handles.sel_par6,'Value',1);
            else
            	set(handles.sel_par6,'Value',0);
            end;
			set(handles.sel_par6,'Enable','on');
			set(handles.sel_par6,'TooltipString',defline);
			set(handles.par6_edit,'String',default);
			set(handles.par6_edit,'Enable','on');
			set(handles.par6_edit,'TooltipString',defline);
        end;
    end;
    % Check for parameter 7
    t1=findstr(tline,'par(7)');
    if ~isempty(t1),
        if tline(1)=='%',
            defline=tline(t1(1)+6:length(tline));
            [mnemo,defline]=strtok(defline);
            handles.model_pars_labels{7}=mnemo;
            [default,defline]=strtok(defline);
            handles.model_defaults(7)=str2double(default);
            [lobound,defline]=strtok(defline);
            handles.model_lower_bounds(7)=str2double(lobound);
            [hibound,defline]=strtok(defline);
            handles.model_upper_bounds(7)=str2double(hibound);
			set(handles.sel_par7,'String',mnemo);
            if enabled>=7,
            	set(handles.sel_par7,'Value',1);
            else
            	set(handles.sel_par7,'Value',0);
            end;
			set(handles.sel_par7,'Enable','on');
			set(handles.sel_par7,'TooltipString',defline);
			set(handles.par7_edit,'String',default);
			set(handles.par7_edit,'Enable','on');
			set(handles.par7_edit,'TooltipString',defline);
        end;
    end;
    % Check for parameter 8
    t1=findstr(tline,'par(8)');
    if ~isempty(t1),
        if tline(1)=='%',
            defline=tline(t1(1)+6:length(tline));
            [mnemo,defline]=strtok(defline);
            handles.model_pars_labels{8}=mnemo;
            [default,defline]=strtok(defline);
            handles.model_defaults(8)=str2double(default);
            [lobound,defline]=strtok(defline);
            handles.model_lower_bounds(8)=str2double(lobound);
            [hibound,defline]=strtok(defline);
            handles.model_upper_bounds(8)=str2double(hibound);
			set(handles.sel_par8,'String',mnemo);
            if enabled>=8,
            	set(handles.sel_par8,'Value',1);
            else
            	set(handles.sel_par8,'Value',0);
            end;
			set(handles.sel_par8,'Enable','on');
			set(handles.sel_par8,'TooltipString',defline);
			set(handles.par8_edit,'String',default);
			set(handles.par8_edit,'Enable','on');
			set(handles.par8_edit,'TooltipString',defline);
        end;
    end;
    % Check for model mode definition 
    t1=findstr(tline,'#extended#');
    if ~isempty(t1),
        if tline(1)=='%',
            handles.model_mode=1;
        end;
    end;
end;
fclose(fid);

handles.model_pars=handles.model_defaults;
handback=handles;