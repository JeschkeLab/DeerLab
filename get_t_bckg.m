function t_bckg=get_t_bckg(handles)

base=handles.APT_kernel;
crosstalk=handles.APT_crosstalk;
trcnorm=handles.APT_norm;
ny=handles.APT_ny;
t=handles.APT_t;

[texp,vexp,zt,phi,imo,dt]=pre_process(handles.t_orig,handles.v_orig,handles.phase,handles.imo,handles.zerotime,handles.dt);
% Determine fit range for background

t_bckg=handles.bckg_start;
t_cutoff=handles.cutoff; %-handles.zerotime;
ttemp=texp-t_bckg*ones(size(texp));
[tt,nofitp0]=min(abs(ttemp));
ttemp=texp-t_cutoff*ones(size(texp));
[tt,pcutoff]=min(abs(ttemp));

imagflag=get(handles.imaginary,'Value');
dcmplx=imagflag*handles.cmplx; % Complex data?

% reference deconvolution for variable-time DEER, normalization for
% constant-time DEER
if handles.ctvt,
    z=handles.A_vexp;
    ref=vexp(1,:);
    sc=max(real(ref));
    sig=vexp(2,:);
    vexp=real(sig)./real(ref)+1i*imag(ref)/sc;
else
	vexp=vexp/max(real(vexp));
end;

[tmi,ztpoi]=min(abs(texp));
texp=texp(ztpoi:pcutoff);
vexp=real(vexp((ztpoi:pcutoff)));
nexp=length(texp);

% Adaptive background correction
nfa=round(0.1*nexp);
if nfa<1, nfa=1; end;
nfe=round(0.6*nexp);
if nfe<5, nfe=5; end;
merit=zeros(1,nfe-nfa);
%     aptmat=zeros(nfe-nfa,149);
for nofitp0=nfa:nfe,
%         disp(sprintf('%s%i','Testing background start at: ',nofitp0));
rt=texp(nofitp0); % time value for start of baseline region
	t_fit=texp(nofitp0:length(texp)); % time window of baseline region
	td_fit=vexp(nofitp0:length(vexp)); % experimental data in this window
	
	% Background fit    
	td_poly=fit_bckg(handles,texp,t_fit,td_fit);
    
	cfac=td_poly(1);
	
	td_exp2=vexp-td_poly; % subtract background
	td_exp2=td_exp2./td_poly; % divide by background, eqn [13]
	%dipevo=td_exp2/cfac2;
	dipevo=td_exp2/max(td_exp2); % normalize
    [m,n]=size(base); % size of kernel
	nt=length(t); % length of time axis
	spc2=zeros(1,m); % initialize distribution
    td=zeros(1,n);
    td(1:length(texp))=dipevo;
	tdx=td.*t; % eqn [21]
	for k=1:m, % sum in eqn [21]
      spc2(k)=spc2(k)+sum(base(k,:).*tdx)/trcnorm(k); 
	end;
	spc3=crosstalk\spc2'; % crosstalk correction, eqn [22]
    merit(nofitp0-nfa+1)=sum(abs(spc3(1:3)));
    tact=texp(nofitp0);
    pstr=sprintf('%s%d%s%6.4f','Optimizing background fit range, Start: ',tact,' ns, Figure of merit: ',abs(spc3(1)));
    set(handles.status_line,'String',pstr);
    drawnow;
end;
[minme,meind]=min(merit);
nofitp1=meind+nfa-1;
t_bckg=texp(nofitp1);
