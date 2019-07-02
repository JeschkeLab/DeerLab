function handback=update_DA(handles)
%
% Update all changed plots in DeerAnalysis
%
% Modified by H.C. Hyde, 2010 (Last updated: 2011-12-20)
%  *Added display of RMS, model type, and distance constraints to workspace
%  *Added call to new function 'stats_analysis' to get F-test results

DEERNet_flag=get(handles.select_deernet,'Value');

bmode = get_bckg_mode(handles);
if strcmp(bmode,'d')
    set(handles.pushbutton_validate_model,'Enable','on');
else
    set(handles.pushbutton_validate_model,'Enable','off');
end

if length(handles.A_low)==length(handles.A_distr)...
        && length(handles.A_high)==length(handles.A_distr)
    cstate = get(handles.error_estimate,'Value');
    if cstate
        if max([max(abs(handles.A_distr-handles.A_low)) max(abs(handles.A_high-handles.A_distr))]) < 10*eps
            set(handles.status_line,'String','Distance distribution error estimate not available (deactivated).');
            set(handles.error_estimate,'Value',0);
        end
    end        
else
    set(handles.error_estimate,'Value',0);
end

see_uncertainty = get(handles.error_estimate,'Value');
see_uncertainty_bckg = get(handles.checkbox_deernet_error_bckg,'Value');
if ~see_uncertainty
    set(handles.error_estimate,'ForegroundColor','r');
    set(handles.error_estimate,'FontWeight','bold');
else
    set(handles.error_estimate,'ForegroundColor','k');
    set(handles.error_estimate,'FontWeight','normal');
end

r = [];
distr = [];

if handles.new_distr
    handles.A_prev_distr = handles.A_curr_distr;
    handles.A_prev_r = handles.A_curr_r;
    handles.A_prev_mode = handles.A_curr_mode;
    handles.new_distr = 0;
end

if handles.new_bckg
    handles.A_prev_bckg_mode = handles.A_curr_bckg_mode;
    handles.A_prev_bckg_details = handles.A_curr_bckg_details;
    handles.A_prev_bckg = handles.A_curr_bckg;
    handles.A_prev_bckg_t = handles.A_curr_bckg_t;    
end

set(handles.text_comparative_bckg,'String',handles.A_prev_bckg_details);
set(handles.text_comparative_bckg,'ForegroundColor',[0.2,0.4,0]);

stats = get(handles.checkbox_statistics,'Value');

full_depth=get(handles.checkbox_model_fit_depth,'Value');

ghost_suppression=get(handles.checkbox_ghost,'Value');

locked=get(handles.checkbox_no_analysis,'Value');
locked=locked*handles.locked_loaded;

smooth_scaled = get(handles.checkbox_smooth_scaled,'Value');

set(handles.constraints,'String','n.a.');

guidance=get(handles.checkbox_guidance,'Value');

if isempty(handles.v_orig), return; end

if handles.dt<=handles.min_dt
    set(handles.dt_minus,'Enable','off');
else
    set(handles.dt_minus,'Enable','on');
end
if handles.dt>=handles.max_dt
    set(handles.dt_plus,'Enable','off');
else
    set(handles.dt_plus,'Enable','on');
end
set(handles.dt_display,'String',sprintf('%i',handles.dt));

% man_bckg_flag=get(handles.manual_bckg,'Value');

dual_flag=get(handles.dual_display,'Value'); % Dual display?
if isempty(handles.B_vexp); dual_flag=0; end
residual_flag=get(handles.residual,'Value'); % show residuum
error_flag=get(handles.error_estimate,'Value'); % show confidence interval for Tikhonov regularization
comp_bckg = get(handles.checkbox_comparative_bckg,'Value');

print=handles.print_request;
handles.print_request=0;

if ~logical(locked)
    [texp,vexp,~,~,imo]=pre_process(handles.t_orig,handles.v_orig,handles.phase,handles.imo,handles.zerotime,handles.dt);
else
    texp=handles.t_orig;
    vexp=handles.v_orig;
    imo=0;
end
imo=max(real(vexp))*imo/max(max(real(handles.v_orig)));

% Determine fit range for background
texp=texp/1000;
t_bckg=handles.bckg_start;
t_cutoff=handles.cutoff;
ttemp=texp-t_bckg*ones(size(texp))/1000;
[~,nofitp0]=min(abs(ttemp));
ttemp=texp-t_cutoff*ones(size(texp))/1000;
[~,pcutoff]=min(abs(ttemp));

imagflag=get(handles.imaginary,'Value');
dcmplx=imagflag*handles.cmplx; % Complex data?

% reference deconvolution for variable-time DEER, normalization for
% constant-time DEER
if handles.ctvt
    ref=vexp(1,:);
    sc=max(real(ref));
    sig=vexp(2,:);
    vexp=real(sig)./real(ref)+1i*imag(ref)/sc;
elseif ~logical(locked) && ~smooth_scaled
	vexp=vexp/max(real(vexp));
end

% Prepare plot window
if print
    figure(1); clf;
    set(gca,'FontSize',24);
else
	axes(handles.original_data);
	cla;
    set(gca,'FontSize',8);
end

handles.texp=texp;
handles.vexp=vexp;

% Plot data sets
hvexp = plot(texp,real(vexp),'k');
set(gca,'FontSize',8);
hold on;
if dual_flag
    B_vexp=real(handles.B_vplot);
    if imagflag
        [~,~,B_vexp,~]=compare_modnorm(texp,real(vexp),handles.B_tplot,real(B_vexp));
        cla;
        set(gca,'FontSize',8);
        plot(texp,real(vexp),'k');
    end
    plot(handles.B_tplot,B_vexp,'b');
end
if logical(dcmplx) && ~dual_flag
	plot(texp,0*vexp,'k:','LineWidth',1.5);
	plot(texp,imag(vexp),'m');
    plot(texp,imo*ones(size(vexp)),'m:');
	ma=max([max(real(vexp)) max(imag(vexp)) max(imo)]);
	mi=min([min(real(vexp)) min(imag(vexp)) min(imo)]);
else
    ma=max(real(vexp));
    mi=min(real(vexp));
    if dual_flag
        if ma<max(B_vexp)
            ma=max(B_vexp); 
        end
        if mi>min(B_vexp)
            mi=min(B_vexp); 
        end
    end
end
sc=ma-mi;
minv=mi-0.1*sc;
maxv=ma+0.1*sc;
axis([min(texp),max(texp),minv,maxv]);
plot([0,0],[minv,maxv],'g');
plot([t_bckg/1000,t_bckg/1000],[minv,maxv],'c');
plot([t_cutoff/1000,t_cutoff/1000],[minv,maxv],'Color',[0.75 0.25 0]);
xlabel('t (µs)');

if comp_bckg
    plot(handles.A_prev_bckg_t,handles.A_prev_bckg,'Color',[0.2,0.4,0],'Linewidth',1);
end

bmode = get_bckg_mode(handles);
if ~strcmp(bmode,'d')
    handles.bckg_request_d = 0;
end

if isfield(handles,'bckg_request_h') && handles.bckg_request_h
    set(handles.bckg_homogeneous,'Value',1);
    set(handles.bckg_none,'Value',0);
elseif isfield(handles,'bckg_request_p') && handles.bckg_request_p
    set(handles.bckg_poly,'Value',1);
    set(handles.bckg_none,'Value',0);
elseif isfield(handles,'bckg_request_e') && handles.bckg_request_e
    set(handles.bckg_exp,'Value',1);
    set(handles.bckg_none,'Value',0);
elseif isfield(handles,'bckg_request_d') && handles.bckg_request_d
    set(handles.bckg_none,'Value',0);
    set(handles.bckg_homogeneous,'Value',0);        
    set(handles.bckg_poly,'Value',0);
    set(handles.bckg_exp,'Value',0);
end

if ~handles.validation_mode && ~handles.bckg_request_d
    % Background fit
    hflag=get(handles.bckg_homogeneous,'Value');
    pflag=get(handles.bckg_poly,'Value');
    eflag=get(handles.bckg_exp,'Value');
    set(handles.bckg_density,'String','n.a.');
    bckg=zeros(size(vexp)); % default background, no correction
    if hflag || pflag || eflag 
        t_fit=texp(nofitp0:pcutoff); % time window of baseline region
        td_fit=real(vexp(nofitp0:pcutoff)); % experimental data in this window
        if min(td_fit)<0.02
            set(handles.bckg_none,'Value',1);
            set(handles.bckg_homogeneous,'Value',0);
            set(handles.bckg_poly,'Value',0);
            set(handles.bckg_exp,'Value',0);
            handles.bckg_fct=zeros(size(texp));
            handles.bckg_rms_value=0;
            set(handles.status_line,'String','Data decay to less than 2% of initial amplitude. Background correction switched off.');
            set(handles.title_original,'ForegroundColor','r');
            handles.bckg_request_h=hflag;
            handles.bckg_request_p=pflag;
            handles.bckg_request_e=eflag;
            hflag=false;
            pflag=false;
            eflag=false;
        else
            if logical(locked)
                bckg=handles.loaded_bckg;
            else
                [bckg,handles]=fit_bckg(handles,texp,t_fit,td_fit);
            end
            diff0=real(vexp(nofitp0:pcutoff))-bckg(nofitp0:pcutoff);
            rms=sqrt(sum(diff0.*diff0)/(length(diff0)-1));
            pstr=sprintf('%8.6f',rms);
            set(handles.bckg_rms,'String',pstr);
            plot(t_fit,bckg(nofitp0:pcutoff),'r','LineWidth',1.5);
            handles.bckg_fct=bckg;
            handles.bckg_rms_value=rms;
            plot(texp,bckg,'r:','LineWidth',1.5);
            handles.bckg_fct = bckg;
            set(handles.title_original,'ForegroundColor','k');
        end
        plot(texp,real(vexp),'k');
    end
    
    if ghost_suppression
        if hflag || pflag || eflag
            raw_ff=real(vexp)./bckg;
        else
            raw_ff=real(vexp);
        end
        cluster=raw_ff.^(1/(handles.spins_per_object-1));
        if hflag || pflag || eflag
            dipevo=cluster-ones(size(cluster));
        else
            dipevo=cluster;
        end

        dipevo=dipevo/max(dipevo);
        cluster=cluster/max(cluster);
    else
        dipevo=real(vexp)-bckg; % subtract background, pure dipolar evolution of local part
        cluster=real(vexp); % without subtraction, retains modulation depth

        if hflag || pflag || eflag
            dipevo=dipevo./bckg; % divide by background, eqn [13]
            cluster=cluster./bckg;
        end
        cluster=cluster/max(cluster);
    end
    [~,ztpoi]=min(abs(texp));

    ttemp=texp-t_cutoff*ones(size(texp))/1000;
    [~,pcutoff]=min(abs(ttemp));

    scale_dipevo=1;
    scale_cluster=1;
    scale_flag=get(handles.renormalize,'Value');
    if scale_flag
        test=dipevo/max(dipevo)-0.5*ones(size(dipevo));
        test=test(ztpoi:length(test));
        flag=0;
        km=round(length(test)/2);
        for k=1:length(test)
            if ~flag
                if test(k)<0
                    km=k+2;
                    flag=1;
                end
            end
        end
        test=test(1:km);
%         disp(km);
%         figure(15); clf;
%         plot(test);

        [~,poi]=min(abs(test));

        apoi=ztpoi-poi;
        epoi=ztpoi+poi;
        
        if apoi<1, apoi=1; end
        if epoi>length(dipevo), epoi=length(dipevo); end
        fit_range=dipevo(apoi:epoi)/max(dipevo);
        time_range=texp(apoi:epoi);
        poly = polyfit(time_range,fit_range,5);
        smoothed=polyval(poly,time_range);
        scale_dipevo=1/max(smoothed);

%         figure(13); clf;
%         plot(time_range,fit_range,'k');
%         hold on;
%         plot(time_range,smoothed,'r');
%         title('Form factor');

        fit_range=cluster(apoi:epoi)/max(cluster);
        time_range=texp(apoi:epoi);
        poly=polyfit(time_range,fit_range,5);
        smoothed=polyval(poly,time_range);
        scale_cluster=1/max(smoothed);

%         figure(14); clf;
%         plot(time_range,fit_range,'k');
%         hold on;
%         plot(time_range,smoothed,'r');
%         title('Cluster');
    end
    
   
    tdip=texp(ztpoi:pcutoff);
    dipevo=dipevo(ztpoi:pcutoff);
    cluster=cluster(ztpoi:pcutoff);
    dipevo_err=bckg(ztpoi:pcutoff);
    bckg=bckg(ztpoi:pcutoff);
    ndip=length(dipevo);
    handles=update_kernel(handles,ndip);

    % Long-pass filtering, if selected
    flag=get(handles.long_pass,'Value');
    if flag
        dipevo=long_pass_filter(handles,tdip,dipevo);
        cluster=long_pass_filter(handles,tdip,cluster);
    end

    dipevo=scale_dipevo*dipevo/max(dipevo);
    cluster=scale_cluster*cluster/max(cluster);
    
    if logical(locked)
        cluster=handles.loaded_ff;
        tdip=handles.loaded_tdip;
    end

    handles.A_dipevo=dipevo;
    handles.A_dipevo_err=ones(size(dipevo_err))./dipevo_err;
    handles.A_cluster=cluster;
    handles.A_bckg=bckg;
    handles.A_tdip=tdip;
    handles.A_depth=1-bckg(1);
elseif handles.validation_mode
    dipevo=handles.A_dipevo;
    cluster=handles.A_cluster;
    bckg=handles.A_bckg;
    tdip=handles.A_tdip;
    dens=handles.bckg_dens;
    plot(tdip,bckg,'r','LineWidth',1.5);
    pstr=sprintf('%6.3f',dens*handles.calib_density);
    if dens>=0
        set(handles.bckg_density,'String',pstr);
    else
        set(handles.bckg_density,'String','n.a.');
    end
    handles.density_value=dens*handles.calib_density;
elseif handles.bckg_request_d
    if ~isempty(handles.A_deernet_bckg)
        cluster = handles.A_deernet_vexp./handles.A_deernet_bckg;
        cluster = cluster/max(cluster);
        tdip = handles.A_deernet_t;
        plen = round(length(cluster)/10);
        p = polyfit(tdip(1:plen),cluster(1:plen),5);
        fit = polyval(p,tdip(1:plen));
        cluster = cluster/max(fit);
        my_offset = 1-handles.A_depth;
        if ghost_suppression
            cluster=cluster.^(1/(handles.spins_per_object-1));
            my_offset = my_offset^(1/(handles.spins_per_object-1)); 
        end
        dipevo = cluster - my_offset;
        dipevo = dipevo/max(dipevo);
        plen = round(length(dipevo)/10);
        p = polyfit(tdip(1:plen),dipevo(1:plen),5);
        fit = polyval(p,tdip(1:plen));
        dipevo = dipevo/max(fit);

        if see_uncertainty_bckg
            delete(hvexp);
            [me,~] = size(handles.deernet_ensemble_bckg);
            for km = 1:me
                plot(handles.A_deernet_t,handles.deernet_ensemble_bckg(km,:),'Color',[1,0.7,0.7],'LineWidth',1);
            end
            plot(texp,real(vexp),'k');
        end
        plot(handles.A_deernet_t,handles.A_deernet_bckg,'r','LineWidth',1.5);
        pstr=sprintf('%6.3f',handles.bckg_dens*handles.calib_density);
        if handles.bckg_dens>=0
            set(handles.bckg_density,'String',pstr);
        else
            set(handles.bckg_density,'String','n.a.');
        end
        handles.A_dipevo = dipevo;
        handles.A_cluster = cluster;
        handles.A_bckg = handles.A_deernet_bckg;
        handles.A_tdip = tdip;
    end
end
pstr=sprintf('%5.3f',handles.A_depth);
set(handles.mod_depth_display,'String',pstr);

if handles.new_bckg
    [bmode,details] = get_bckg_mode(handles);
    handles.A_curr_bckg_mode = bmode;
    handles.A_curr_bckg_details = details;
    handles.A_curr_bckg = handles.A_bckg;
    handles.A_curr_bckg_t = handles.A_tdip;
end

handles.new_bckg = false;

if exist('dipevo','var')
    test=dipevo-0.5*ones(size(dipevo));
    flag=0;
    km=round(length(test)/2);
    for k=1:length(test)
        if ~flag
            if test(k)<0
                km=k+2;
                flag=1;
            end
        end
    end
    test=test(1:km);
    [~,poi]=min(abs(test));
    tcross=tdip(poi);
    dist_est=5*(tcross/0.48)^(1/3);
    % disp(sprintf('%s%4.2f%s','Half value: ',tcross,' µs'));
    % disp(sprintf('%s%4.2f%s','Distance estimate: ',dist_est,' nm'));
    across=cluster(poi);
    handles.dist_est=dist_est;

    if print
        figure(2); clf;
        set(gca,'FontSize',24);
    else
        axes(handles.dipolar_evolution);
        cla;
        set(gca,'FontSize',8);
    end

    zfm=handles.zf-1;
    dipevo2=[private_hamming(dipevo,1) zeros(1,zfm*length(dipevo))];
    dipevo2(1)=dipevo2(1)/2;
    dipspc=real(fftshift(fft(dipevo2)));
    fmin=-length(dipevo)/(2*max(tdip));          % s. Schweiger/Jeschke book p. 106
    fmax=(length(dipevo)-1)/(2*max(tdip));
    faxges=linspace(fmin,fmax,length(dipspc));
    handles.A_spc=dipspc;
    handles.A_ny=faxges;

    flag=get(handles.dip_time_domain,'Value');
    if flag
        if handles.bckg_request_d && exist('cluster','var') && handles.select_deernet.Value
            plot(tdip,real(cluster),'k');
            set(gca,'FontSize',8);
            scd = max(real(cluster))-min(real(cluster));
            mind=min(real(cluster))-0.1*scd;
            maxd=max(real(cluster))+0.1*scd;
            hold on;
            if handles.updated
                plot(handles.A_deernet_t,handles.A_deernet_ff,'r');
                scd = max(handles.A_deernet_ff)-min(handles.A_deernet_ff);
                mind=min(handles.A_deernet_ff)-0.1*scd;
                maxd=max(handles.A_deernet_ff)+0.1*scd;
                if residual_flag
                    cla;
                    residual = real(cluster)-handles.A_deernet_ff;
                    plot(handles.A_deernet_t,residual,'k');
                    scd = max(residual)-min(residual);
                    mind=min(residual)-0.1*scd;
                    maxd=max(residual)+0.1*scd;
                end
            end
            axis([min(handles.A_deernet_t),max(handles.A_deernet_t)/handles.zoom,mind,maxd]);
        else
            handles.text_form_factor.String = 'Form factor';
            handles.text_form_factor.ForegroundColor = [0,0,0];
            plot(tdip,real(cluster),'k');
            set(gca,'FontSize',8);
            hold on;
            plot(tcross,across,'bo','MarkerFaceColor','b');
            text(tcross+0.02*max(tdip),across,sprintf('%4.1f%s',dist_est,' nm'),'Color','b');
            if dual_flag
                B_cluster=handles.B_cluster;
                if imagflag
                    [~,~,B_cluster]=compare_modnorm_linear(tdip,cluster,handles.B_tdip,real(B_cluster),handles.A_depth,handles.B_depth);
                end
                plot(handles.B_tdip,B_cluster,'b');
            end
            xlabel('t (µs)');
            mac=max(cluster);
            mic=min(cluster);
            if dual_flag
                if max(B_cluster)>mac
                    mac=max(B_cluster);
                end
                if min(B_cluster)<mic
                    mic=min(B_cluster);
                end
            end
            scd=mac-mic;
            mind=mic-0.1*scd;
            maxd=mac+0.1*scd;
            axis([min(tdip),max(tdip)/handles.zoom,mind,maxd]);
        end
    else
        mis=min(dipspc);
        mas=max(dipspc);
        if dual_flag
            B_spc=handles.B_spc;
            if imagflag
                B_spc=B_spc*max(dipspc)/max(B_spc);
            end
            if min(B_spc)<mis
                mis=min(B_spc);
            end
            if max(B_spc)>mas
                mas=max(B_spc);
            end
            plot(handles.B_ny,handles.B_spc,'b');
            hold on;
        end
        plot(faxges,dipspc,'k');
        xlabel('f (MHz)');
        scd=mas-mis;
        mind=mis-0.1*scd;
        maxd=mas+0.1*scd;
        axis([fmin/handles.zoom,fmax/handles.zoom,mind,maxd]);
    end
    hold on

    r = handles.A_r;
    distr = handles.A_distr;
    sim = handles.A_sim;


    APT_flag=get(handles.select_APT,'Value');
    if APT_flag && ~logical(locked)
        [r,distr,sim]=APT(handles,tdip,dipevo);
        if ~handles.updated, handles.mask=ones(size(r)); end
        handles.APT=distr;
        handles.r_APT=r;
        handles.A_r=r;
        handles.A_distr=distr;
        handles.mask=ones(size(distr));
        handles.updated=1;
        handles.APT_sim=sim;
        handles.A_sim=sim;
    end

    if logical(locked)
        handles.APT=distr;
        handles.r_APT=r;
        sim=handles.loaded_sim;
        distr=handles.loaded_distr;
        r=handles.loaded_rax;
        handles.A_r=r;
        handles.A_distr=distr;
        handles.A_sim=sim;
        handles.updated=1;
        handles.mask=ones(size(distr));
    end
    
    % New code for log display and fit
    % theoretically nice, but numerically problematic
    % if handles.log_Tikh,
    %     logsim=log(sim);
    %     logcluster=log(handles.A_cluster);
    %     logdipevo=logcluster;
    %     sc=sum(logdipevo.*logdipevo)/sum(logdipevo.*logsim);
    %     logsim=sc*logsim;
    %     figure(13); clf;
    %     plot(tdip,logdipevo,'k');
    %     hold on;
    %     plot(tdip,logsim,'r');
    % 	axes(handles.dipolar_evolution);
    % end;

    drmin=min(r);
    drmax=max(r);
    pstr=sprintf('%s%3.1f%s%3.1f%s','Range: (',drmin,',',drmax,') nm');
    set(handles.distr_range,'String',pstr);

    sim0=sim-0.99*ones(size(sim)); % standard modulation depth 1 %
    
    if DEERNet_flag && handles.updated
        deernet_dipevo = handles.A_deernet_ff - handles.A_deernet_bckg(1);
        sim0 = 0.01*deernet_dipevo/max(deernet_dipevo);
    end

    spcflag=get(handles.dip_spectrum,'Value');
    if handles.updated % && ~handles.bckg_request_d
        if ~full_depth
            modsim=ones(size(sim))-sim;
            modexp=ones(size(cluster))-cluster;
            sc=sum(modexp.*modexp)/sum(modsim.*modexp);
            sim=ones(size(modsim))-sc*modsim;
        end
        difference=sim-cluster;
        rms=sqrt(sum(difference.*difference)/(length(difference)-1));
        handles.fit_rms_value=rms;
        if rms < handles.best_rmsd
            handles.best_rmsd = rms;
        end
        pstr=sprintf('%8.6f',rms);
        set(handles.distr_rms,'String',pstr);
        if ~handles.validation_mode
            % modulation depth lambda is 0.01*sc
            if sc>99
                numspin=-1;
                numspinstr='n.a.';
            else
                numspin=1+log((1-0.01*sc/handles.moddepth_suppression))/handles.calib_numspins;
                numspinstr=sprintf('%4.2f',numspin);
            end
            handles.n_spins=numspin;
            set(handles.num_spins,'String',numspinstr);
        else
            numspin = handles.n_spins;
        end
        sim0(1)=sim0(1)/2;
        simspc=real(fftshift(fft([private_hamming(sim0,1) zeros(1,zfm*length(sim0))])));
        sim0(1)=2*sim0(1);
        handles.A_simspc=simspc;
        if spcflag
            sc0=sum(dipspc.*dipspc)/sum(simspc.*dipspc);
            simspc=sc0*simspc;
            scs=max(simspc)-min(simspc);
            mins=min(simspc)-0.1*scs;
            maxs=max(simspc)+0.1*scs;
            mind=min([mind mins]);
            maxd=max([maxd maxs]);
            plot(faxges,simspc,'r','LineWidth',1.5);
            plot(faxges,dipspc,'k');
            if residual_flag
                cla;
                set(gca,'FontSize',8);
                residual=dipspc-simspc;
                plot(faxges,residual,'k');
                scs=max(residual)-min(residual);
                mind=min(residual)-0.1*scs;
                maxd=max(residual)+0.1*scs;
            end
            axis([fmin/handles.zoom,fmax/handles.zoom,mind,maxd]);
        elseif ~handles.select_deernet.Value
            if sum(handles.mask)<length(handles.mask)
                plot(tdip,handles.mask_sim,'g','LineWidth',1.5);
                difference=handles.mask_sim-cluster;
                rms=sqrt(sum(difference.*difference))/(length(difference)-1);
                handles.fit_rms_value=rms;
                rmsstr=sprintf('%s%8.6f','r.m.s. deviation with supression: ',rms);
                set(handles.status_line,'String',rmsstr);
            end
            plot(tdip,sim,'r','LineWidth',1.5);
            plot(tdip,cluster,'k');

            plot(tcross,across,'bo','MarkerFaceColor','b');
            text(tcross+0.02*max(tdip),across,sprintf('%4.1f%s',dist_est,' nm'),'Color','b');
            if residual_flag
                cla;
                set(gca,'FontSize',8);
                residual=cluster-sim;
                sresidual=residual;
                if sum(handles.mask)<length(handles.mask)
                    sresidual=cluster-handles.mask_sim;
                    plot(tdip,sresidual,'g','LineWidth',1.5);
                end
                plot(tdip,residual,'k');
                mas=max([max(residual) max(sresidual)]);
                mis=min([min(residual) min(sresidual)]);
                scd=mas-mis;
                mind=mis-0.1*scd;
                maxd=mas+0.1*scd;
                axis([0,max(tdip)/handles.zoom,mind,maxd]);
            end
        end
        % Determine suggested cutoff 
        tr=tdip(6:length(tdip)-5);
        mvavg=tr;
        for k=1:length(mvavg)
            mvavg(k)=norm(difference(k:k+10));
        end
        [min_noise,poi]=min(mvavg);
        cut=length(mvavg);
        while mvavg(cut-1)>6*min_noise && cut>poi
            cut=cut-1;
        end
        tcut=tdip(cut+10);
        handles.cutoff_suggestion=1000*tcut;
        if ~spcflag
            plot([tcut,tcut],[min(cluster),max(cluster)],'Color',[0.75 0.25 0]);
        end
    end
end

if stats
    % New code (HCH): Display RMS, model type, and distance constraints in workspace 
    % --------
    if handles.updated && (handles.Tikh_updated || handles.model_updated)

        % STATISTICAL ANALYSIS
        % The RMS threshold based on the F-test for the current fit is displayed in
        %   the workspace. It is ONLY valid if the current model is unconstrained and
        %   the global minimum was effectively reached!
        model_stats = stats_analysis(handles,rms);
        handles.model_stats = model_stats;
        % model_stats

        % Display model type and parameter summary in workspace and if applicable, 2-gaussian fit distance ratio
        if handles.model_updated
            if (strcmp(handles.user_model,'Two_Gaussians') || strcmp(handles.user_model,'Two_Gaussians_hom') || strcmp(handles.user_model,'Two_Rice3d')) && handles.confit
                k1 = handles.model_pars_kcon(1);
                k2 = handles.model_pars_kcon(2);
                if handles.nconstraints > 0
                    % Display a 2-gaussian or 2-rice model with nonlinear parameter constraints
                    fprintf(1,'\n%s (%s %s) [%s%d,%s%d,%s%d] %s%s%0.4f\n','Fit:',handles.user_model,'+NLcon','p=',model_stats.nfp,...
                                 'v=',model_stats.ndf,'n=',model_stats.nobs,handles.nlconstraint_label,'=',handles.model_pars(k2)/handles.model_pars(k1));
                else
                    % Display a standard 2-gaussian or 2-rice model (unconstrained)
                    fprintf(1,'\n%s (%s) [%s%d,%s%d,%s%d] %s%s%0.4f\n','Fit:',handles.user_model,'p=',model_stats.nfp,...
                                 'v=',model_stats.ndf,'n=',model_stats.nobs,handles.nlconstraint_label,'=',handles.model_pars(k2)/handles.model_pars(k1));
                end
            elseif ~isempty(model_stats)    
                % Display all other standard models
                fprintf(1,'\n%s (%s) [%s%d,%s%d,%s%d]\n','Fit:',handles.user_model,'p=',model_stats.nfp,'v=',model_stats.ndf,'n=',model_stats.nobs);  
            end
        end

        % Display fit and statistical results in workspace. Notation: finv(_,_,_) = finv(1-P,ndf_num,ndf_denom)
        % USER MODEL-specific
        if handles.model_updated && ~isempty(model_stats)
           fprintf(1,'%s %9.7f\n','RMS  =',model_stats.rms);
            if handles.confit && handles.nconstraints == 0
                model_disp = 1;   %workspace display version {1,2}: show F distribution notation vs. Matlab command for FINV function
                switch model_disp
                  case 1
                    fprintf(1,'\n%s %9.7f [%s(%s%d,%s%d,%s%0.3f)]\n',...
                        'RMSt =',model_stats.rms_thresh,'by F-test for F','v1=',model_stats.Fndf_numer,'v2=',model_stats.Fndf_denom,'1-',model_stats.P);
                  case 2
                   fprintf(1,'\n%s %9.7f [%s(%s%0.3f,%d,%d)]\n',...
                        'RMSt =',model_stats.rms_thresh,'by F-test for finv','1-',model_stats.P,model_stats.Fndf_numer,model_stats.Fndf_denom);
                end
            end
        end

        % TIKHONOV-specific
        if handles.Tikh_updated && ~isempty(model_stats)
            % Display Tikhonov regularization ("eff" designates effective values!)
            fprintf(1,'\n%s (%s%s) [%s%d,%s%d,%s%d]\n','Fit:','Tikhonov: Reg.par.=',num2str(handles.regpar),'p_eff=',model_stats.nfp,'v_eff=',model_stats.ndf,'n=',model_stats.nobs);
            fprintf(1,'%s %9.7f\n','RMS  =',model_stats.rms);
        end

        % FOR ALL FITS: AICc (display the corrected Akaike Information Criterion score for this fit)
        if (handles.model_updated || handles.Tikh_updated) && ~isempty(model_stats)
            fprintf(1,'\n%s %9.7f\n','AICc =',model_stats.AICc);
        end

    end
    % --------
end

fit_flag=get(handles.select_model,'Value');

if exist('dipevo','var') 
    if fit_flag && ~handles.updated
        [r,distr]=APT(handles,tdip,dipevo);
        handles.APT=distr;
        handles.r_APT=r;
        if length(handles.model_sim)~=length(cluster)
            handles=sim_user_model(handles);
        end
        model_sim=handles.model_sim;
        model_sim0=model_sim-0.99*ones(size(model_sim));
        modsim=ones(size(model_sim))-model_sim;
        modexp=ones(size(cluster))-cluster;
        sc=sum(modexp.*modexp)/sum(modsim.*modexp);
        model_sim=ones(size(modsim))-sc*modsim;
        model_sim0(1)=model_sim0(1)/2;
        simspc=real(fftshift(fft([private_hamming(model_sim0,1) zeros(1,3*length(sim0))])));
        handles.A_simspc=simspc;
        if spcflag
            sc0=sum(dipspc.*dipspc)/sum(simspc.*dipspc);
            simspc=sc0*simspc;
            scs=max(simspc)-min(simspc);
            mins=min(simspc)-0.1*scs;
            maxs=max(simspc)+0.1*scs;
            mind=min([mind mins]);
            maxd=max([maxd maxs]);
            plot(faxges,simspc,'r:','LineWidth',1.5);
            plot(faxges,dipspc,'k');
            axis([fmin/handles.zoom,fmax/handles.zoom,mind,maxd]);
        else
            plot(tdip,model_sim,'r:','LineWidth',1.5);
            plot(tdip,cluster,'k');
            difference=model_sim-cluster;
            rms=sqrt(sum(difference.*difference)/(length(difference)-1));
            handles.fit_rms_value=rms;
            if rms < handles.best_rmsd
                handles.best_rmsd = rms;
            end
            pstr=sprintf('%8.6f',rms);
            set(handles.distr_rms,'String',pstr);
        end
    end
end

if print
    figure(3); clf;
    set(gca,'FontSize',24);
else
	axes(handles.distance_distribution);
    set(handles.distance_distribution,'NextPlot','replace');
	cla;
    set(gca,'FontSize',8);
end

LC_flag=get(handles.L_curve,'Value');

if ~isempty(r) && ~isempty(distr) && handles.updated && ~LC_flag
	rmin=handles.rmin;
	rmax=handles.rmax;
    mind=1;
    maxd=10;
    eflag=get(handles.distr_expand,'Value');
    if eflag
        mind=rmin;
        maxd=rmax;
    end
    dist_intg=sum(distr); % unscaled integral of distance distribution
    if handles.n_spins>0
        sc_dist=real((handles.n_spins-1)/dist_intg); % scaling factor for distance distribution 
    else
        sc_dist=1;
    end
    sc=max(distr)-min(distr);
	minv=min(distr)-0.1*sc;
	maxv=max(distr)+0.1*sc;
    if error_flag && length(handles.A_low)==length(distr) && length(handles.A_high)==length(distr)
        minv = min(handles.A_low)-0.1*sc;
        maxv = max(handles.A_high)+0.1*sc;
    end
    maxt=max(tdip);
    sc=(maxt/2)^(1/3);
    rel_shape=sc*3;
    rel_mean_width=sc*4;
    rel_mean=sc*5;
    recognize=sc*6;
    if guidance
        rectangle('Position',[min(r),sc_dist*minv,rel_shape-min(r),sc_dist*(maxv-minv)],'EdgeColor','none','FaceColor',[0.9,1,0.9]);
        set(gca,'FontSize',8);
        rectangle('Position',[rel_shape,sc_dist*minv,rel_mean_width-rel_shape,sc_dist*(maxv-minv)],'EdgeColor','none','FaceColor',[1,1,0.8]);
        rectangle('Position',[rel_mean_width,sc_dist*minv,rel_mean-rel_mean_width,sc_dist*(maxv-minv)],'EdgeColor','none','FaceColor',[1,0.9,0.7]);
        rectangle('Position',[rel_mean,sc_dist*minv,recognize-rel_mean,sc_dist*(maxv-minv)],'EdgeColor','none','FaceColor',[1,0.8,0.8]);
    end   
    sc=max(distr)-min(distr);
    if sum(handles.mask) ~= length(handles.mask)
        plot(r,sc_dist*handles.mask.*distr,'g','LineWidth',1.5);
        set(gca,'FontSize',8);
    end
	hold on;
    
    dlow=handles.A_low;
    dhigh=handles.A_high;
    if error_flag && length(dlow)==length(distr) && length(dhigh)==length(distr)
%         for k=1:length(r)
%             plot([r(k) r(k)],sc_dist*[dlow(k) dhigh(k)],'Color',[0.65 0.65 0.65],'Linewidth',0.5);
%         end;
        [mr,nr] = size(r);
        if mr > nr
            r = r';
        end
        hunc = fill([r fliplr(r)],sc_dist*[dlow fliplr(dhigh)],[0.6,0.6,0.6]);
        set(hunc,'EdgeColor','none','FaceAlpha',0.75);
        plot(r,sc_dist*handles.mean_distr,'k','LineWidth',1.5);
        sc=max(dhigh)-min(dlow);
        maxv=max(dhigh)+0.1*sc;
        minv=min(dlow)-0.1*sc;
    end

    if dual_flag
        B_distr=handles.B_distr;
        if ~imagflag
            B_distr=B_distr*max(distr)/max(B_distr);
        end
        plot(handles.B_r,sc_dist*B_distr,'b','LineWidth',1.5);
        if minv>min(B_distr)-0.1*sc
            minv=min(B_distr)-0.1*sc;
        end
        if maxv<max(B_distr)+0.1*sc
            maxv=max(B_distr)+0.1*sc;
        end
    end
    if handles.theor
        th_distr=handles.th_distr;
        th_distr=max(sc_dist*distr)*th_distr/max(th_distr);
        plot(handles.th_r,th_distr,'c','LineWidth',1.5);
    end
    if handles.checkbox_comparative.Value && ~isempty(handles.A_prev_r)
        scp = max(sc_dist*distr)/max(handles.A_prev_distr);
        plot(handles.A_prev_r,scp*handles.A_prev_distr,'Color',[0.2,0.4,0],'LineWidth',1);
        set(handles.text_comparative,'String',handles.A_prev_mode);
        set(handles.text_comparative,'ForegroundColor',[0.2,0.4,0]);
    else
        set(handles.text_comparative,'String',' ');
        set(handles.text_comparative,'ForegroundColor','k');
    end
    plot(r,sc_dist*distr,'k','LineWidth',1.5);
    handles.A_curr_r = r;
    handles.A_curr_distr = distr;
    sel_APT = get(handles.select_APT,'Value');
    if sel_APT
        handles.A_curr_mode = 'APT';
    end
    sel_Tikh = get(handles.select_Tikhonov,'Value');
    if sel_Tikh
        handles.A_curr_mode = sprintf('Tikhonov %14.5f',handles.regpar);
    end
    sel_deernet = get(handles.select_deernet,'Value');
    if sel_deernet
        net_set = get(handles.popupmenu_deernet,'Value');
        net_sets = get(handles.popupmenu_deernet,'String');
        net_set = net_sets{net_set};
        handles.A_curr_mode = sprintf('DEERNet %s',net_set);
    end
    sel_model = get(handles.select_model,'Value');
    if sel_model
        my_model = get(handles.user_model_list,'Value');
        my_models = get(handles.user_model_list,'String');
        my_model = strtrim(my_models(my_model,:));
        handles.A_curr_mode = sprintf('Model %s',my_model);
    end
    
    xlabel('{\itr} (nm)');
    ylabel('{\itP} (nm^{-1})');
	
	plot([rmin,rmin],sc_dist*[minv,maxv],'b');
	plot([rmax,rmax],sc_dist*[minv,maxv],'m');
    if recognize<maxd
        maxd=recognize;
    end
	axis([mind,maxd,sc_dist*minv,sc_dist*maxv]);
else
    sc_dist=1;
end

if fit_flag && ~handles.updated && ~LC_flag
    model_distr=handles.model_distr;
    plot(handles.r_APT,handles.APT,'k:','LineWidth',1);
    set(gca,'FontSize',8);
    hold on;
    model_distr=model_distr*max(handles.APT)/max(model_distr);
    plot(handles.model_r,model_distr,'r:','LineWidth',1.5);
    xlabel('{\itr} (nm)');
    ylabel('{\itP} (nm^{-1})');
	
	rmin=handles.rmin;
	rmax=handles.rmax;
    minv=min(handles.APT);
    maxv=max(handles.APT);
    scd=maxv-minv;
    minv=minv-0.1*scd;
    maxv=maxv+0.1*scd;
    mind=1;
    maxd=10;
    eflag=get(handles.distr_expand,'Value');
    if eflag
        mind=rmin;
        maxd=rmax;
    end
	plot([rmin,rmin],[minv,maxv],'b');
	plot([rmax,rmax],[minv,maxv],'m');
	axis([mind,maxd,minv,maxv]);
end

if LC_flag
    h = plot(handles.Lcurve_rho,handles.Lcurve_eta,'k.');
    set(h,'PickableParts','none');
    set(gca,'FontSize',8);
    xlabel(gca,'log(||S-KP||) (misfit)');
    ylabel(gca,'log(||LP||) (roughness)');
    hold on;
    poi=handles.regpar_sel;
    if handles.regpar_sel == handles.regpar_opt_Lc
        mark = [1,0,0];
    elseif handles.regpar_sel == handles.regpar_opt_AIC
        mark = [1,0,0];
    elseif handles.regpar_sel == handles.regpar_opt_GCV
        mark = [1,0,0];
    else
        mark = [0,0.5,1];
    end
    h = plot(handles.Lcurve_rho(poi),handles.Lcurve_eta(poi),'o','Color',mark);
    set(h,'PickableParts','none');
    set(gca,'ButtonDownFcn',@(hObject,eventdata)DeerAnalysis('LCurve_ButtonDownFcn',hObject,eventdata,guidata(hObject)));
else
    set(gca,'ButtonDownFcn',@(hObject,eventdata)DeerAnalysis('empty_ButtonDownFcn',hObject,eventdata,guidata(hObject)));    
end

if handles.updated
	rmin=handles.rmin;
	rmax=handles.rmax;
	if rmin<min(r), rmin=min(r); end
	if rmax>max(r), rmax=max(r); end
	rnum=round((rmax-rmin)/0.02);
	rcut=linspace(rmin,rmax,rnum);
	distcut=interp1(r,distr,rcut,'pchip',0);
	mom=moment_analysis_vec(rcut,distcut);
	handles.moments=mom;
	handles.r_MA=mom(1);
	handles.sr_MA=mom(2);
    pstr=sprintf('%5.2f',mom(1));
    set(handles.r_mean,'String',pstr);
    pstr=sprintf('%5.2f',real(sqrt(mom(2))));
    set(handles.sr,'String',pstr);
    % determine and display constraints
    h=6.626176e-34;
    g=2.0059;
    muB=9.274078e-24;
    mu0=4*pi*1e-7;
    intg_distr=cumsum(distcut);
%     figure(13); clf;
%     plot(rcut,intg_distr,'k');
%     hold on;
%     plot([rcut(1),rcut(length(rcut))],[0.1*sum(distcut),0.1*sum(distcut)],'c');
%     plot([rcut(1),rcut(length(rcut))],[0.9*sum(distcut),0.1*sum(distcut)],'c');
    [~,min_poi]=min(abs(intg_distr-0.1*sum(distcut)*ones(size(intg_distr))));
    handles.min_constraint=rcut(min_poi);
    [~,max_poi]=min(abs(intg_distr-0.9*sum(distcut)*ones(size(intg_distr))));
    handles.max_constraint=rcut(max_poi);
    % test if upper constraint makes sense, (see Eq. (2.8) in: 
    % G. Jeschke, Ye. Polyhach, Phys. Chem. Chem. Phys. 9, 1895-1910 (2007)
    ttest=1e6*8*pi*h*handles.max_constraint^3*1e-27/(g^2*muB^2*mu0);
    pstr=sprintf('%s%4.2f%s','[',handles.min_constraint,', ');
    if max(tdip)<0.5*ttest
        handles.max_constraint=1000;
        pstr=sprintf('%s%s',pstr,'n.a.]');
        set(handles.manual_bckg,'ForegroundColor','k');
    else
        set(handles.manual_bckg,'ForegroundColor','r');
        pstr=sprintf('%s%4.2f%s',pstr,handles.max_constraint,']');
    end
    set(handles.constraints,'String',pstr);
    if ~LC_flag
        plot([handles.min_constraint,handles.min_constraint],sc_dist*[minv,maxv],'c:');
        plot([handles.max_constraint,handles.max_constraint],sc_dist*[minv,maxv],'c:');
    end
else
    set(handles.r_mean,'String','n.a.');
    set(handles.sr,'String','n.a.');
    set(handles.constraints,'String','n.a.');
end

handback=handles;
% Update handles structure
guidata(gcbo, handles);