function [mode,details] = get_bckg_mode(handles)

mode = '';
details = '';
nm = get(handles.bckg_none,'Value');
if nm
    mode = 'n';
    details = 'none';
else
    hm = get(handles.bckg_homogeneous,'Value');
    if hm
        mode = 'h';
        details = 'hom';
        homdimstr = get(handles.bckg_dim_edit,'String');
        homdim = str2double(homdimstr);
        details = sprintf('%s(%4.2f)',details,homdim);
    else
        pm = get(handles.bckg_poly,'Value');
        if pm
            mode = 'p';
            details = 'poly';
            polystr = get(handles.bckg_poly_order,'String');
            polyorder = str2double(polystr);
            details = sprintf('%s(%i)',details,polyorder);
        else
            em = get(handles.bckg_exp,'Value');
            if em
                mode = 'e';
                details = 'experimental';
            else
                dm = get(handles.radiobutton_deernet_bckg,'Value');
                if dm
                    mode = 'd';
                    details = handles.comp_net_set;
                end
            end
        end
    end
end