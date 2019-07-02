function handback=process_list(handles)
%
% Load data set in the selected format and save it in data structure handles
%

autosave=get(handles.autosave,'Value');

handback=handles;

reset_flag=get(handles.reset_on_load,'Value');
set(handles.reset_on_load,'Value',0);

Lflag=get(handles.select_L_curve,'Value');
set(handles.select_L_curve,'Value',0);

data_format=get(handles.format_winepr,'Value'); % 0 for Elexsys or ASCII, 1 for WINEPR
flag=get(handles.format_ascii,'Value');
if flag, data_format=2; end;

set(handles.main_figure,'Pointer','watch');

phasing=get(handles.autophase,'Value');

[fname,pname]=uigetfile('*.txt;*.lst','Open list of filenames');
if isequal(fname,0)|isequal(pname,0), 
    return; 
end;

dots=findstr('.',fname);
ext_pos=length(fname);
if length(dots)>0, ext_pos=dots(length(dots))-1; end;
list_name=fname(1:ext_pos); % name without extension


fid = fopen( [ pname fname], 'r');% par
[S,sn]= fscanf(fid,'%c');
fclose(fid);

cd(pname);
currdir=pwd;

%Parse string S, the delimiter is ascii 10, this is the token delimiter
D = setstr(10);
VB=[];
while(length(S) > 0 )
[token,S] = strtok(S,D);
VB = str2mat(VB,token);
end

nf=0;

[m,n]=size(VB);

% Find out, how many filenames in list
for k=1:m,
    act=deblank(VB(k,:));
    if sum(size(act)) & act(1)~='#'>0,
        nf=nf+1;
    end;
end;

if nf==0, return; end;

list_m1=zeros(1,nf);
list_m2=zeros(1,nf);
list_m3=zeros(1,nf);
weight=zeros(1,nf);
rstd=linspace(1.0,10.0,2000);
mean_dipevo=zeros(1,2048);
arr_vexp=zeros(nf,2048);
arr_texp=zeros(nf,2048);
max_nexp=0;
max_poi_arr=0;
min_poi_arr=16384;
cutpoi=0;

mean_dist=zeros(1,length(rstd));

success=0; % flag that is set to 1 if at least one spectrum is analyzed
poi=0; % counter for spectra that can actually be analyzed by the selected method

for k=1:m,
    act=deblank(VB(k,:));
    if sum(size(act)) & act(1)~='#'>0,
        set(handles.data_set_A,'String',act);
        bas_name=act;
        handles.source_file=[pname act];
		switch data_format,
            case 0, % Elexsys
                [x,y,z,vb]=get_elexsys(act);
            case 1, % WINEPR
                [x,y,z,vb]=get_WINEPR(act);
            case 2, % ASCII
                dset=load_ascii(act,2);
                [mtest,ntest]=size(dset); % determine size of loaded data array
                x=dset(:,handles.ascii_t);
                z=dset(:,handles.ascii_real);
                if ntest>=handles.ascii_imag,
                    z=z+i*dset(:,handles.ascii_imag);
                end;
                x=real(x')-i*imag(x');
                z=z';
                y=[];
                vb=[];
		end;
		[mtest,ntest]=size(z); % determine size of loaded data array
		if mtest==2,
            handles.ctvt=1;
            dtype='Variable-time DEER. ';
		else,
            handles.ctvt=0;
            dtype='Constant-time DEER. ';
		end;
		cflag=round(1000*sum(abs(imag(z)))/sum(abs(real(z)))); % Determine, if imaginary part is significant
		if cflag,
            handles.cmplx=1;
            dform=' complex ';
            handles.dcmplx=1;
		else,
            handles.cmplx=0;
            handles.dcmplx=0;
            dform=' real ';
		end;
		nexp=length(x);
		pstr=sprintf('%d',nexp);
		set(handles.status_line,'String',[dtype pstr dform 'data points.']);
		handles.bas_name=bas_name;
		currdir=pwd;
		handles.A_vb=vb;
		handles.maxt=max(x);
		phi=handles.phase;
        handles.dt=x(2)-x(1);

		handles.A_texp=x;
		handles.A_vexp=z;
		handles.A_vb=vb;

        % Determine fit range for background
		texp=(handles.A_texp-handles.zerotime*ones(size(handles.A_texp)))/1000;
		t_bckg=handles.bckg_start;
		t_cutoff=handles.cutoff-handles.zerotime;
		ttemp=texp-handles.zerotime*ones(size(texp))/1000;
		[tt,ztpoi]=min(abs(ttemp));
		ttemp=texp-t_bckg*ones(size(texp))/1000;
		[tt,nofitp0]=min(abs(ttemp));
		ttemp=texp-t_cutoff*ones(size(texp))/1000;
		[tt,pcutoff]=min(abs(ttemp));
		
		dcmplx=handles.cmplx; % Complex data?
        
		phaseflag=get(handles.autophase,'Value');
		if phaseflag,
            handles=get_phase(handles);
		end;

		% phase correction, for variable-time DEER also reference deconvolution
		phi=handles.phase;
		if handles.ctvt,
            z=handles.A_vexp;
            ref=z(1,:)*exp(i*phi);
            sc=max(real(ref));
            sig=z(2,:)*exp(i*phi);
            vexp=real(sig)./real(ref)+i*imag(ref)/sc;
		else,
			vexp=handles.A_vexp*exp(i*phi);
			vexp=vexp/max(real(vexp));
		end;
        
        % Store experimental spectrum between zero time and cutoff time
        vexp_cut=vexp(ztpoi:pcutoff);
        texp_cut=texp(ztpoi:pcutoff);
        if length(vexp_cut)>max_poi_arr, max_poi_arr=length(vexp_cut); end;
        if length(vexp_cut)<min_poi_arr, min_poi_arr=length(vexp_cut); end;
        cutpoi=cutpoi+1;
        arr_vexp(cutpoi,1:length(vexp_cut))=vexp_cut;
        arr_texp(cutpoi,1:length(vexp_cut))=texp_cut;
        
		axes(handles.original_data);
		cla;
		% Background fit
		hflag=get(handles.bckg_homogeneous,'Value');
		pflag=get(handles.bckg_poly,'Value');
		eflag=get(handles.bckg_exp,'Value');
		set(handles.bckg_density,'String','n.a.');
		bckg=zeros(size(vexp)); % default background, no correction
		if hflag | pflag | eflag, 
			t_fit=texp(nofitp0:pcutoff); % time window of baseline region
			td_fit=real(vexp(nofitp0:pcutoff)); % experimental data in this window
			
			[bckg,handles]=fit_bckg(handles,texp,t_fit,td_fit);
			diff0=real(vexp(nofitp0:pcutoff))-bckg(nofitp0:pcutoff);
			rms=sqrt(sum(diff0.*diff0))/(length(diff0)-1);
			pstr=sprintf('%8.6f',rms);
			set(handles.bckg_rms,'String',pstr);
            plot(t_fit,bckg(nofitp0:pcutoff),'r','LineWidth',1.5);
            handles.bckg_fct=bckg;
            handles.bckg_rms_value=rms;
            plot(texp,bckg,'r:','LineWidth',1.5);
		end;
		plot(texp,real(vexp),'k');
        
        sc=max(real(vexp))-min(real(vexp));
        axis([min(texp),max(texp),min(real(vexp))-0.1*sc,max(real(vexp))+0.1*sc]);

        
		dipevo=real(vexp)-bckg; % subtract background, pure dipolar evolution of local part
		cluster=real(vexp); % without subtraction, retains modulation depth
		
		if hflag | pflag | eflag,
			cfac=bckg(1);
            cfac2=0;
            if cfac>0,
				cfac2=(1-cfac)/cfac;
            end;
			
			dipevo=dipevo./bckg; % divide by background, eqn [13]
			cluster=cluster./bckg;
		end;
		cluster=cluster/max(cluster);
		
		[tmi,ztpoi]=min(abs(texp));
		
		ttemp=texp-t_cutoff*ones(size(texp))/1000;
		[tt,pcutoff]=min(abs(ttemp));
		
		tdip=texp(ztpoi:pcutoff);
		dipevo=dipevo(ztpoi:pcutoff);
		cluster=cluster(ztpoi:pcutoff);
		ndip=length(dipevo);
		handles=update_kernel(handles,ndip);
		
		% Long-pass filtering, if selected
		flag=get(handles.long_pass,'Value');
		if flag,
            dipevo=long_pass_filter(handles,tdip,dipevo);
            cluster=long_pass_filter(handles,tdip,cluster);
		end;
		
		dipevo=dipevo/max(dipevo);
		cluster=cluster/max(cluster);
		
		handles.A_dipevo=dipevo;
		handles.A_cluster=cluster;
		handles.A_tdip=tdip;

		dipevo2=[private_hamming(dipevo,1) zeros(1,3*length(dipevo))];
		dipspc=real(fftshift(fft(dipevo2)));
		dt=tdip(2)-tdip(1);
		fmin=-(1+1/length(dipspc))/(2*dt);          % s. Schweiger/Jeschke book p. 106
		fmax=1/(2*dt);
		faxges=linspace(fmin,fmax,length(dipspc));
		handles.A_spc=dipspc;
		handles.A_ny=faxges;
        
        if length(dipevo)>max_nexp, max_nexp=length(dipevo); end;
        
        APT_flag=get(handles.select_APT,'Value');
		if APT_flag,
            [r,distr,sim]=APT(handles,tdip,dipevo);
            if ~handles.updated, handles.mask=ones(size(r)); end;
            handles.A_r=r;
            handles.A_distr=distr;
            handles.APT_sim=sim;
            handles.A_sim=sim;
		end;

        Tikh_flag=get(handles.select_Tikhonov,'Value');
		if Tikh_flag,
            handles=fit_Tikhonov_new(handles);
		end;

        model_flag=get(handles.select_model,'Value');
		if model_flag,
			handles=fit_user_model(handles);
		end;

        poi=poi+1;
        rms=1.0e6;
        sim=handles.A_sim;
        cluster=handles.A_cluster;
		axes(handles.dipolar_evolution);
		cla;
		if length(handles.A_sim)>0,
            sim0=sim-0.99*ones(size(sim)); % standard modulation depth 1 %
            simspc=real(fftshift(fft([private_hamming(sim0,1) zeros(1,3*length(sim0))])));
            handles.A_simspc=simspc;
			modsim=ones(size(sim))-sim;
			modexp=ones(size(cluster))-cluster;
			sc=sum(modexp.*modexp)/sum(modsim.*modexp);
			sim=ones(size(modsim))-sc*modsim;
			numspin=1+log((1-0.01*sc))/handles.calib_numspins;
			handles.n_spins=numspin;
			numspinstr=sprintf('%4.2f',numspin);
			set(handles.num_spins,'String',numspinstr);
            diff=sim-cluster;
            rms=sqrt(sum(diff.*diff))/(length(diff)-1);
            handles.fit_rms_value=rms;
            pstr=sprintf('%8.6f',rms);
			set(handles.distr_rms,'String',pstr);
            success=1;
            plot(tdip,sim,'r','LineWidth',1.5);
		end;
        plot(tdip,cluster,'k');
        sc=max(cluster)-min(cluster);
        axis([min(tdip),max(tdip),min(cluster)-0.1*sc,max(cluster)+0.1*sc]);

        if autosave,
            handles=save_result(handles,bas_name,pname);
        end;
        

        r=handles.A_r;
        distr=handles.A_distr;

		axes(handles.distance_distribution);
		cla;
        plot(r,distr,'k','LineWidth',1.5);
	    xlabel('r (nm)');
        sc=max(distr)-min(distr);
        axis([min(r),max(r),min(distr-0.1*sc),max(distr)+0.1*sc]);
        
        rmin=handles.rmin;
        rmax=handles.rmax;
		rnum=round((rmax-rmin)/0.02);
        
		rcut=linspace(rmin,rmax,rnum);
		distcut=interp1(r,distr,rcut,'pchip',0);
		
		% Moment analysis
		mom=moment_analysis_vec(rcut,distcut);
		handles.moments=mom;
		handles.rmean=mom(1);
		handles.sigr=sqrt(mom(2));        


        stddist=get_std_distr(r,distr,rstd);
        vari=rms^2;
        weight(poi)=1/vari;
        mean_dist=mean_dist+stddist/vari;
        mean_dipevo(1:length(dipevo))=mean_dipevo(1:length(dipevo))+dipevo/vari;
		rmin=handles.rmin;
		rmax=handles.rmax;
        if rmin<min(r), rmin=min(r); end;
        if rmax>max(r), rmax=max(r); end;
		rnum=round((rmax-rmin)/0.02);

        rcut=linspace(rmin,rmax,rnum);
		distcut=interp1(r,handles.A_distr,rcut,'pchip',0);
		
		% Moment analysis
		mom=moment_analysis_vec(rcut,distcut);
        list_m1(poi)=mom(1);
        list_m2(poi)=mom(2);
        list_m3(poi)=mom(3);          
    end;
end;

arr_vexp=real(arr_vexp(1:cutpoi,1:max_poi_arr));

% Make data array with all data sets normalized to the modulation depth of
% the first data set

arr_modnorm=arr_vexp(:,1:min_poi_arr);

for k=2:cutpoi,
    [tdipn,trace0,trace1,stddev]=compare_modnorm_old(arr_texp(1,:),arr_vexp(1,:),arr_texp(k,:),arr_vexp(k,:));
    arr_modnorm(k,:)=trace1(1:min_poi_arr);
end;
arr_modnorm(1,:)=trace0(1:min_poi_arr);
arr_modnorm=[texp_cut(1:min_poi_arr)' arr_modnorm'];

% Make comparison matrix

cmpmat=zeros(cutpoi,cutpoi);
warning off MATLAB:polyfit:RepeatedPointsOrRescale
for k=1:cutpoi,
    for l=1:cutpoi,
        [tstd,trace0,trace1,stddev]=compare_modnorm_old(arr_texp(l,:),arr_vexp(l,:),arr_texp(k,:),arr_vexp(k,:));
        if k~=l,
            xax=linspace(1,length(trace0),length(trace0));
            [p0,s]=polyfit(xax,trace0,13);
            decay0=polyval(p0,xax);
            [p1,s]=polyfit(xax,trace1,13);
            decay1=polyval(p1,xax);
            ctrace0=trace0-decay0;
            ctrace1=trace1-decay1;
            diff1=trace0-trace1;
            rms1=sqrt(sum(diff1.*diff1));
            diff2=ctrace0-ctrace1;
            rms2=sqrt(sum(diff2.*diff2));
            if rms2>0,
                cmpmat(k,l)=rms1/rms2-1;
            else,
                cmpmat(k,l)=0;
            end;
        else,
            cmpmat(k,l)=0;
        end;
    end;
end;

weight=weight(1:poi);
list_m1=list_m1(1:poi);
list_m2=list_m2(1:poi);
list_m3=list_m3(1:poi);
mean_dist=mean_dist/sum(weight);
mean_dipevo=mean_dipevo/sum(weight);
mean_dipevo=mean_dipevo(1:max_nexp);
dt=texp(2)-texp(1);
tstd=linspace(0,(length(mean_dipevo)-1)*dt,length(mean_dipevo));

meanr=sum(list_m1.*weight)/sum(weight);
stdr=sqrt(var(list_m1,weight));
mean_sig=sum(sqrt(list_m2).*weight)/sum(weight);
stdstd=sqrt(var(sqrt(list_m2),weight));
mean_m3=sum(list_m3.*weight)/sum(weight);
stdm3=sqrt(var(list_m3,weight));

save_list_result(handles,list_name,rstd,mean_dist,tstd,mean_dipevo,meanr,stdr,mean_sig,stdstd,mean_m3,stdm3,arr_modnorm,cmpmat);

% Update handles structure
guidata(handles.main_figure, handles);

set(handles.main_figure,'Pointer','arrow');

set(handles.reset_on_load,'Value',reset_flag);
set(handles.select_L_curve,'Value',Lflag);

handback=handles;
