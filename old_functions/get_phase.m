function handback=get_phase(handles)
% Determine phase correction

texp=(handles.t_orig-handles.zerotime*ones(size(handles.t_orig)))/1000;
t_bckg=handles.bckg_start;
t_cutoff=handles.cutoff-handles.zerotime;
ttemp=texp-t_bckg*ones(size(texp))/1000;
[tt,nofitp0]=min(abs(ttemp));
ttemp=texp-t_cutoff*ones(size(texp))/1000;
[tt,pcutoff]=min(abs(ttemp));

if handles.ctvt,
    tr=handles.v_orig(1,:);
else
    tr=handles.v_orig;
end;
tr=tr(nofitp0:pcutoff);
last=tr(end);
phi0=atan2(imag(last),real(last));
im_off_0=0;
v=[phi0,im_off_0];
pstart=round(length(tr)/8); % use only last 7/8 of data for phase/offset correction
v=v(1);
v=fminsearch(@rms_phi,v,[],tr(pstart:end));
phi=v(1);
if sum(real(tr*exp(1i*phi)))<0, phi=phi+pi; end;
handles.phase=phi;
handles.imo=0;
pstr=sprintf('%6.1f',phi*180/pi);
set(handles.phase_edit,'String',pstr);

% Update handles structure
guidata(handles.main_figure, handles);
handback=handles;
