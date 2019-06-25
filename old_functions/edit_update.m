function [v,handback]=edit_update(handles,hObject,vmin,vmax,vdef,fstr,intgflag)
% Update an edit field
% with check whether value is in proper range (vmin,vmax) and reset to
% default value vdef if necessary, format string fstr for sprintf
% determines appearance
% if intgflag=1, value is an integer and is rounded
% the corrected value is the output parameter
%

set(handles.status_line,'String','Value edited.');
set(handles.status_line,'ForegroundColor','b');
v=str2double(get(hObject,'String'));
if intgflag, v=round(v); end;
% Protect against wrong inputs
if isnan(v),
    v=vdef;
    set(handles.status_line,'String','### Not a number. Resetting to default value. ###');
end;
if v<vmin, 
    v=vmin;
    set(handles.status_line,'String','### Smaller than lower limit. Resetting to default value. ###');
end;
if v>vmax, 
    v=vmax;
    set(handles.status_line,'String','### Larger than upper limit. Resetting to default value. ###');
end;
pstr=sprintf(fstr,v);
set(hObject,'String',pstr);
% Update handles structure
guidata(hObject, handles);
handback=handles;
