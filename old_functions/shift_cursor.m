function [v,handback]=shift_cursor(handles,hObject,vmin,vmax,fstr,v0,increment),
% Shifts a cursor by increment
% with check whether value stays in proper range (vmin,vmax) format string fstr for sprintf
% determines appearance in associated edit field (hObject)
% the new value is the output parameter
%

set(handles.status_line,'String','Value incremented/decremented.');
set(handles.status_line,'ForegroundColor','b');
v=v0+increment;
% Protect against wrong inputs
if isnan(v),
    v=v0;
    set(handles.status_line,'String','### Not a number. Resetting to default value. ###');
end;
if v<vmin, 
    v=v0;
    set(handles.status_line,'String','### Smaller than lower limit. Resetting to default value. ###');
end;
if v>vmax, 
    v=v0;
    set(handles.status_line,'String','### Larger than upper limit. Resetting to default value. ###');
end;
pstr=sprintf(fstr,v);
set(hObject,'String',pstr);
% Update handles structure
guidata(hObject, handles);
handback=handles;
