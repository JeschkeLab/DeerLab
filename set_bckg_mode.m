function set_bckg_mode(handles,bmode)

switch bmode
    case 'n'
        set(handles.bckg_none,'Value',1);
        set(handles.bckg_homogeneous,'Value',0);
        set(handles.bckg_poly,'Value',0);
        set(handles.bckg_exp,'Value',0);
        set(handles.radiobutton_deernet_bckg,'Value',0);
    case 'h'
        set(handles.bckg_none,'Value',0);
        set(handles.bckg_homogeneous,'Value',1);
        set(handles.bckg_poly,'Value',0);
        set(handles.bckg_exp,'Value',0);
        set(handles.radiobutton_deernet_bckg,'Value',0);
    case 'p'
        set(handles.bckg_none,'Value',0);
        set(handles.bckg_homogeneous,'Value',0);
        set(handles.bckg_poly,'Value',1);
        set(handles.bckg_exp,'Value',0);
        set(handles.radiobutton_deernet_bckg,'Value',0);
    case 'e'
        set(handles.bckg_none,'Value',0);
        set(handles.bckg_homogeneous,'Value',0);
        set(handles.bckg_poly,'Value',0);
        set(handles.bckg_exp,'Value',1);
        set(handles.radiobutton_deernet_bckg,'Value',0);
    case 'd'
        set(handles.bckg_none,'Value',0);
        set(handles.bckg_homogeneous,'Value',0);
        set(handles.bckg_poly,'Value',0);
        set(handles.bckg_exp,'Value',0);
        set(handles.radiobutton_deernet_bckg,'Value',1);
end
guidata(handles.bckg_none,handles);