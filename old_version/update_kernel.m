function handback=update_kernel(handles,ndip),
%
% Checks, if correct APT kernel is loaded
% if not, the kernel is updated
%

pot=0; rem=ndip; % find minimum power of 2 which is still larger than number of data points 
while rem>1,
    pot=pot+1;
    rem=rem/2;
end;
ksize=2^pot; % actual kernel size
if ksize ~= handles.kernel_size,
	p_string=sprintf('%d',ksize); % display
	fname=sprintf('%s%s','kernel',p_string); % generate filename
	load(fname);
	handles.APT_kernel=base;
	handles.APT_crosstalk=crosstalk;
	handles.APT_norm=tnorm;
	handles.APT_ny=ny;
	handles.APT_t=t;
    handles.kernel_size=ksize;
end;
handback=handles;

% Update handles structure
guidata(handles.main_figure, handles);
