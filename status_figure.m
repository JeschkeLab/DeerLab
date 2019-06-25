function active = status_figure(fraction,modal,id,title)
% function active = status_figure(fraction,modal,id)
%
% Display a status window that informs the user on progress of a
% computation and allows the user to gently cancel the computation by
% closing this window
%
% Input:
%
% fraction  fraction (0...1) of the computation that is completed, 0, or 
%           empty input reset figure to a newly
%           started computation, figure is closed for fraction = 1
% modal     optional flag defining if the figure is modal, preventing the
%           user from accessing other parts of a complex GUI, defaults to
%           false
% id        optional identifier of the figure for applications that use
%           several status windows, defaults to 1
%
% Ouptut:
%
% active    flag that informs the user if the figure was close, true if the
%           computation should remain active, false if it should be
%           cancelled
%
% ### status_figure defines a persistent variable ###
%
% status_parameters     array of structures with elements corresponding to
%                       distinct status figures, structure has fields
%   .start_time         time when the computation was started
%   .fig_handle         handle of the status figure
%   .text_percent       handle of the percentage string
%   .text_timeleft      handle of the string for the time left
%   .active             flag that indicates whether this computation is
%                       active

persistent status_parameters

if ischar(fraction),
    title = fraction;
end;

if ~exist('fraction','var') || isempty(fraction) || ischar(fraction),
    fraction = 0;
end;

if fraction > 1 || fraction < 0,
    error('status_figure:argChk','Fraction outside allowed range.');
end;

if ~exist('modal','var') || isempty(modal),
    modal = false;
end;

if ~exist('id','var') || isempty(id),
    id = 1;
end;

if ~exist('title','var'),
    title = 'Close to interrupt computation';
end;

active = true;

if ~exist('status_parameters','var') || isempty(status_parameters),
    status_parameters(id).start_time = clock;
    status_parameters(id).active = 1;
    status_parameters(id).fig_handle = [];
    status_parameters(id).text_percent = [];
    status_parameters(id).text_timeleft = [];
    fraction = 0;
end;

if ~isempty(status_parameters(id).fig_handle),
    if ~ishandle(status_parameters(id).fig_handle),
        status_parameters(id).active = 0;
        status_parameters(id).fig_handle = [];
        active = false;
        return
    end;
end;

if fraction < eps,
    status_parameters(id).start_time = clock;
end;

% create figure if required
if (fraction == 0) && isempty(status_parameters(id).fig_handle),
    s=[320 80];
    t=get(0,'ScreenSize');
    fh=figure('HandleVisibility','on','MenuBar',...
        'none','Name',title,'NumberTitle',...
        'off','Resize','off','Position',[floor((t(3:4)-s)/2) s],'Tag',...
        sprintf('StatusFigure%i',id),'ToolBar','none','Visible','on');
    status_parameters(id).fig_handle = fh;
    t1h = uicontrol(fh,'Style','text',...
                'FontWeight','bold',...
                'HorizontalAlignment','left',...
                'String','Completed:',...
                'Position',[15 40 100 20]);
    t2h = uicontrol(fh,'Style','text',...
                'FontWeight','bold',...
                'HorizontalAlignment','left',...
                'String','Time left:',...
                'Position',[15 20 100 20]);
    t3h = uicontrol(fh,'Style','text',...
                'FontWeight','normal',...
                'HorizontalAlignment','left',...
                'String','Completed:',...
                'Position',[90 40 100 20]);
    t4h = uicontrol(fh,'Style','text',...
                'FontWeight','normal',...
                'HorizontalAlignment','left',...
                'String','Time left:',...
                'Position',[90 20 100 20]);
    status_parameters(id).text_percent = t3h;
    status_parameters(id).time_left = t4h;
    set(t3h,'String',sprintf('%4.2f%%',100*fraction));
    set(t4h,'String',sprintf('%s',timestr(status_parameters(id).start_time,fraction)));
    return
end;

% delete figure if computation has completed
% create figure if required
if (fraction == 1) && ~isempty(status_parameters(id).fig_handle),
    delete(status_parameters(id).fig_handle);
    active = false;
    status_parameters(id).active = 0;
    status_parameters(id).fig_handle = [];
    return
end;

if isempty(status_parameters(id).fig_handle) || ~ishandle(status_parameters(id).fig_handle),
    error('status_figure:no_figure','Status figure not yet initialized.');
end;

set(status_parameters(id).text_percent,'String',sprintf('%4.2f%%',100*fraction));
set(status_parameters(id).time_left,'String',sprintf('%s',timestr(status_parameters(id).start_time,fraction)));

if modal,
    set(status_parameters(id).fig_handle,'WindowStyle','modal');
end;
    
figure(status_parameters(id).fig_handle);

function time_string = timestr(start_time,fraction)

time_now = clock;
time_passed = etime(time_now,start_time);
if fraction < eps,
    time_string = 'n/a';
    return
else
    time_left = (1-fraction)*time_passed/fraction;
    d = floor(time_left/(24*60*60));
    time_left = time_left - d*24*60*60;
    h = floor(time_left/(60*60)); 
    time_left = time_left - h*60*60;
    min = floor(time_left/(60)); 
    time_left = time_left - min*60;
    time_string = sprintf('%i s',round(time_left));
    if min > 0,
        time_string = sprintf('%i min %s',min,time_string);
    end;
    if h > 0,
        time_string = sprintf('%i h %s',h,time_string);
    end;
    if d > 0,
        time_string = sprintf('%i d %s',d,time_string);
    end;
end;
