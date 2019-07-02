function handles=save_result(handles,bas_name,path)
%
% Saves a log file, listing the parameters of data analysis and the result of the moment analysis, 
% the original data and background correction,
% the experimental and simulated dipolar evolution functions,
% and the distance distribution
% 
%
% Content and file name convention:
% - filenames are derived from outname
% - the method string is appended to the name of the data set
% - up to six ASCII output files are generated (depending on available data), 
%   which are distinguished by another
%   appendix to the file name:
%   _res    main parameters of the analysis and moment analysis of the
%           distribution
%   _bckg   background fit, three or four columns, time axis, primary experimental data,
%           background function, if the input data were complex, the fourth column is
%           the imaginary part of the primary data after phase correction
%   _fit    fit of the dipolar evolution function, three columns,
%           time axis, experimental dipolar evolution function, simulated
%           dipolar evolution function
%           from version 2015 on, if range suppression is activated, fourth
%           column contains the fit with range suppression
%   _spc    fit of the dipolar spectrum, three columns,
%           frequency axis, experimental dipolar spectrum, simulated
%           dipolar spectrum
%   _distr  distance distribution, two columns, distance axis and
%           distribution, for Tikhonov reg. 4 columns, minimum and maximum
%           of confidence interval (take confidence interval with a grain
%           of salt)
%           from version 2015 on, if range suppression is activated, fifth
%           column (with Tikhonov regularization) or third column (without)
%           contains the distance distribution with range suppression
%   _Lcurve Tikhonov L curve, three columns, rho, eta, and corresponding
%           regularization parameters
% - files _distr, _bckg, and _res have the extension '.dat', file _res has
%   the extension '.txt'
% - time axes are in microseconds, distance axes are in nanometers
%
% Modified by H.C. Hyde, 2010 
%  *Added constraint text for parameters 6,7,8 in case a user model is used
%   with constraints 'kmin' or 'kmax'. 
%  *Added extra fit output section: includes <r2>/<r1> ratio for 2-gaussian
%   or <nu2>/<nu1> for 2-rice model fits, number of data observations, 
%   fitted parameters, and degrees of freedom (including Tikhonov fit).
%
% last edit: G. Jeschke, 10.3.2015

fname_distr=[path bas_name '_distr.dat'];
fname_bckg=[path bas_name '_bckg.dat'];
fname_fit=[path bas_name '_fit.dat'];
fname_spc=[path bas_name '_spc.dat'];
fname_Lcurve=[path bas_name '_Lcurve.dat'];
fname_res=[path bas_name '_res.txt'];

texp=handles.texp;
nexp=length(texp);
vexp=handles.vexp;
tdip=handles.A_tdip;
dipevo=handles.A_cluster;
sim=handles.A_sim;
ny=handles.A_ny;
spc=handles.A_spc;
simspc=handles.A_simspc;
r=handles.A_r;
distr=handles.A_distr;
dlow=handles.A_low;
dhigh=handles.A_high;
bckg=handles.bckg_fct;
rms_bckg=handles.bckg_rms_value;
rmin=handles.rmin;
rmax=handles.rmax;

DEERNet_flag=get(handles.select_deernet,'Value');

if ~isempty(distr)
    % Moment analysis for range of interest
    % Cutout of distance distribution between cursors
    if rmin<min(r), rmin=min(r); end
    if rmax>max(r), rmax=max(r); end
    rnum=round((rmax-rmin)/0.02);
    rcut=linspace(rmin,rmax,rnum);
    distcut=interp1(r,distr,rcut,'pchip',0);
    mom_roi=moment_analysis_vec(rcut,distcut);
    
    % Moment analysis of full distribution, if different from range of interest
    if abs(rmin-min(r))>0.05 || abs(max(r)-rmax)>0.05
        mom_full=moment_analysis_vec(r,distr);
    end
    
end

% Normalize simulated dipolar evolution function, dipolar spectrum, and distance distribution
if ~isempty(bckg)
    if DEERNet_flag
        bckg = nan(size(texp));
        for k = 1:length(texp)
            [mi,poi] = min(abs(tdip-texp(k)));
            if mi < 2*eps
                bckg(k) = handles.A_deernet_bckg(poi);
            end
        end
    end
	data2=[texp' real(vexp') bckg' -imag(vexp')];
	save(fname_bckg,'data2','-ascii');
end
if ~isempty(sim)
	modsim=ones(size(sim))-sim;
	modexp=ones(size(dipevo))-dipevo;
	sc=sum(modexp.*modexp)/sum(modsim.*modexp);
    sim=ones(size(modsim))-sc*modsim;
    if DEERNet_flag
        sim = handles.A_deernet_ff;
    end
	data3=[tdip' dipevo' sim'];
    if sum(handles.mask)<length(handles.mask)
        data3 = [data3 handles.mask_sim'];
    end
	save(fname_fit,'data3','-ascii');
end
if ~isempty(simspc)
    sc0=sum(spc.*spc)/sum(simspc.*spc);
    simspc=sc0*simspc;
	data4=[ny' spc' simspc'];
	save(fname_spc,'data4','-ascii');
end
if ~isempty(distr)
	sc_dist=(handles.n_spins-1)/sum(distr); % scaling factor for distance distribution 
	distr=sc_dist*distr;
    [md,nd] = size(r);
    if md > nd
        r = r';
    end
    [md,nd] = size(distr);
    if md > nd
        distr = distr';
    end
    data1=[r' distr'];
    [md,nd] = size(dlow);
    if md > nd
        dlow = dlow';
    end
    [md,nd] = size(dhigh);
    if md > nd
        dhigh = dhigh';
    end
    if length(dlow)==length(distr) && length(dhigh)==length(distr)
        dlow=sc_dist*dlow;
        dhigh=sc_dist*dhigh;
        data1=[r' distr' dlow' dhigh'];
    end
    if sum(handles.mask)<length(handles.mask)
        data1 = [data1 (distr.*handles.mask)'];
    end
    save(fname_distr,'data1','-ascii');
end

m=length(handles.regpars);
if m>1 % if L_curve exists
    [ml,nl] = size(handles.Lcurve_rho);
    if ml < nl
        data5=[handles.Lcurve_rho' handles.Lcurve_eta' handles.regpars'];
    else
        data5=[handles.Lcurve_rho handles.Lcurve_eta handles.regpars];
    end
	save(fname_Lcurve,'data5','-ascii');
end

rms=handles.fit_rms_value;

% Determine analysis method
method='Approximate Pake transformation';
APT_flag=get(handles.select_APT,'Value');
Tikh_flag=get(handles.select_Tikhonov,'Value');
if Tikh_flag
    method='Tikhonov regularization';
	Lcurve_flag=get(handles.select_L_curve,'Value');
    if Lcurve_flag
        method=[method ' with L curve computation'];
    end
end
model_flag=get(handles.select_model,'Value');
if model_flag
    method=['Fit by model: ' handles.user_model];
end
if DEERNet_flag
    method = 'DEERNet';
end

% Determine background model
bckg_model='No';
noflag=1;
pflag=get(handles.bckg_poly,'Value');
if pflag, bckg_model='Polynomial'; end
hflag=get(handles.bckg_homogeneous,'Value');
if hflag
    bckg_model='Homogeneous ';
    pstr=num2str(handles.hom_dim);
    bckg_model=[bckg_model pstr ' dimensional'];
    dflag=get(handles.bckg_fit_dim,'Value');
end
uflag=get(handles.bckg_exp,'Value');
if uflag
    bckg_model='Fixed polynomial (experimental)';
end

wfile=fopen(fname_res,'w+');
fprintf(wfile,'%s%s%s\n','>>> DEER analysis of data set: ',bas_name,' <<<');
fprintf(wfile,'%s%s%s\n\n','- ',method,' -');
fprintf(wfile,'%s\n','### Description of data set ###');
fprintf(wfile,'%s\n%s\n','Source file:',handles.source_file);
if handles.cmplx, dtype=' complex'; else, dtype=' real'; end
if handles.ctvt, etype='Variable-time DEER. '; else, etype='Constant-time DEER. '; end
fprintf(wfile,'%s%d%s%s\n',etype,nexp,dtype,' data points');
fprintf(wfile,'\n%s\n','### Pre-processing ###');
fprintf(wfile,'%s%d%s\n','Time shift: ',handles.zerotime,' ns');
fprintf(wfile,'%s%d%s\n','Cutoff at : ',handles.cutoff,' ns');
fprintf(wfile,'%s%5.1f%s\n','Phase correction: ',180*handles.phase/pi,'°');
fprintf(wfile,'%s%s\n',bckg_model,' background correction');
fprintf(wfile,'%s%d%s\n','Starting at : ',handles.bckg_start,' ns');
if pflag || uflag
	fprintf(wfile,'%s\n','Background polynomial (for logarithmic data):');
    fprintf(wfile,'%s','(');
    for k=1:length(handles.polynomial)
        fprintf(wfile,'%0.4g',handles.polynomial(k));
        if k<length(handles.polynomial)
            fprintf(wfile,'%s',', ');
        else
            fprintf(wfile,'%s\n',')');
        end
    end
end
fprintf(wfile,'%s%9.6f\n','r.m.s. error of baseline fit: ',rms_bckg);
fprintf(wfile,'%s%s\n','Background density: ',get(handles.bckg_density,'String'));
if ~noflag
	fprintf(wfile,'%s%d%s\n','Background fitting starts at : ',handles.bckg_start,' ns');
end
lp_flag=get(handles.long_pass,'Value');
if lp_flag
    fprintf(wfile,'%s%5.2f%s\n','Distance long-pass filter at ',handles.longpass_min,' nm');
end
if APT_flag
    fprintf(wfile,'%s%5.2f%s\n','Distance-domain smoothing filter width: ',handles.DDS,' nm');
end
fprintf(wfile,'\nModulation depth:  %5.3f\n',handles.A_depth);
if Tikh_flag || model_flag
    exflag=get(handles.exci_bandwidth_corr,'Value');
    if exflag
        fprintf(wfile,'\n%s%5.2f%s\n','Correction for excitation bandwidth of ',handles.bandwidth,' MHz');
    end
end
if ~isempty(distr)
	fprintf(wfile,'\n%s\n','### Distance distribution ###');
    fprintf(wfile,'%s%s\n','Number of spins: ',get(handles.num_spins,'String'));
	fprintf(wfile,'%s%9.6f\n','r.m.s. error of fit: ',rms);
	fprintf(wfile,'%s%5.2f%s%5.2f%s\n','Range of interest from ',rmin,' to ',rmax,' nm');
	fprintf(wfile,'%s%5.2f%s\n','Mean distance  <r> ',mom_roi(1),' nm');
	fprintf(wfile,'%s%5.2f%s\n','Std. dev. sigma(r) ',real(sqrt(mom_roi(2))),' nm');
	fprintf(wfile,'%s%5.2f%s\n','3rd moment         ',mom_roi(3),' nm');
	if abs(rmin-min(r))>0.05 || abs(max(r)-rmax)>0.05
        fprintf(wfile,'\n%s%5.2f%s%5.2f%s\n','Full range from ',min(r),' to ',max(r),' nm');
		fprintf(wfile,'%s%5.2f%s\n','Mean distance  <r> ',mom_full(1),' nm');
		fprintf(wfile,'%s%5.2f%s\n','Std. dev. sigma(r) ',real(sqrt(mom_full(2))),' nm');
		fprintf(wfile,'%s%5.2f%s\n\n','3rd moment         ',mom_full(3),' nm');
	end
	if Tikh_flag
        fprintf(wfile,'%s%s\n','Regularization parameter: ',num2str(handles.regpar));
	end
    if model_flag
        fprintf(wfile,'%s%s\n','Model parameters for ',handles.user_model);
        enstate=get(handles.sel_par1,'Enable');
        enflag=strncmp(enstate,'on',2);
        if enflag
            fprintf(wfile,'%s%s%s',get(handles.sel_par1,'String'),' = ',get(handles.par1_edit,'String'));
            fitflag=get(handles.sel_par1,'Value');
            if fitflag
                fprintf(wfile,'%s\n',', fitted');
            else
                fprintf(wfile,'%s\n',', fixed');
            end
        end
        enstate=get(handles.sel_par2,'Enable');
        enflag=strncmp(enstate,'on',2);
        if enflag
            fprintf(wfile,'%s%s%s',get(handles.sel_par2,'String'),' = ',get(handles.par2_edit,'String'));
            fitflag=get(handles.sel_par2,'Value');
            if fitflag
                fprintf(wfile,'%s\n',', fitted');
            else
                fprintf(wfile,'%s\n',', fixed');
            end
        end
        enstate=get(handles.sel_par3,'Enable');
        enflag=strncmp(enstate,'on',2);
        if enflag
            fprintf(wfile,'%s%s%s',get(handles.sel_par3,'String'),' = ',get(handles.par3_edit,'String'));
            fitflag=get(handles.sel_par3,'Value');
            if fitflag
                fprintf(wfile,'%s\n',', fitted');
            else
                fprintf(wfile,'%s\n',', fixed');
            end
        end
        enstate=get(handles.sel_par4,'Enable');
        enflag=strncmp(enstate,'on',2);
        if enflag
            fprintf(wfile,'%s%s%s',get(handles.sel_par4,'String'),' = ',get(handles.par4_edit,'String'));
            fitflag=get(handles.sel_par4,'Value');
            if fitflag
                fprintf(wfile,'%s\n',', fitted');
            else
                fprintf(wfile,'%s\n',', fixed');
            end
        end
        enstate=get(handles.sel_par5,'Enable');
        enflag=strncmp(enstate,'on',2);
        if enflag
            fprintf(wfile,'%s%s%s',get(handles.sel_par5,'String'),' = ',get(handles.par5_edit,'String'));
            fitflag=get(handles.sel_par5,'Value');
            if fitflag
                fprintf(wfile,'%s\n',', fitted');
            else
                fprintf(wfile,'%s\n',', fixed');
            end
        end
        stats = get(handles.checkbox_statistics,'Value');
        if stats && handles.fit_constrained
            % ---New code(HCH): Begin to write nonlinar constraint info to saved results file---
            % Assumes that nonlin constraint bounds are located at parameter index >= 6
            % (continues until end of file)
            enstate=get(handles.sel_par6,'Enable');
            enflag=strncmp(enstate,'on',2);
            if enflag
                fprintf(wfile,'%s%s%s',get(handles.sel_par6,'String'),' = ',get(handles.par6_edit,'String'));
                fitflag=get(handles.sel_par6,'Value');
                if any(strcmp(get(handles.sel_par6,'String'),{'kmin','kmax'}))
                    if fitflag
                        fprintf(wfile,'%s\n',', active constraint');
                    else
                        fprintf(wfile,'%s\n',', inactive constraint');
                    end
                else
                    if fitflag
                        fprintf(wfile,'%s\n',', fitted');
                    else
                        fprintf(wfile,'%s\n',', fixed');
                    end
                end
            end
            enstate=get(handles.sel_par7,'Enable');
            enflag=strncmp(enstate,'on',2);
            if enflag
                fprintf(wfile,'%s%s%s',get(handles.sel_par7,'String'),' = ',get(handles.par7_edit,'String'));
                fitflag=get(handles.sel_par7,'Value');
                if any(strcmp(get(handles.sel_par6,'String'),{'kmin','kmax'}))
                    if fitflag
                        fprintf(wfile,'%s\n',', active constraint');
                    else
                        fprintf(wfile,'%s\n',', inactive constraint');
                    end
                else
                    if fitflag
                        fprintf(wfile,'%s\n',', fitted');
                    else
                        fprintf(wfile,'%s\n',', fixed');
                    end
                end
            end
            enstate=get(handles.sel_par8,'Enable');
            enflag=strncmp(enstate,'on',2);
            if enflag
                fprintf(wfile,'%s%s%s',get(handles.sel_par8,'String'),' = ',get(handles.par8_edit,'String'));
                fitflag=get(handles.sel_par8,'Value');
                if any(strcmp(get(handles.sel_par6,'String'),{'kmin','kmax'}))
                    if fitflag
                        fprintf(wfile,'%s\n',', active constraint');
                    else
                        fprintf(wfile,'%s\n',', inactive constraint');
                    end
                else
                    if fitflag
                        fprintf(wfile,'%s\n',', fitted');
                    else
                        fprintf(wfile,'%s\n',', fixed');
                    end
                end
            end
        end % HC Hyde code
    end
end

if Tikh_flag && isfield(handles,'model_stats') && ~isempty(handles.model_stats)
    %For Tikhonov fit, write number of data observations, fitted parameters
    %and degrees of freedom
    fprintf(wfile,'\n%s\n','### Extra fit output ###');
    fprintf(wfile,'%s %0d\n','Data observations      :',handles.model_stats.nobs);
    fprintf(wfile,'%s %0d\n','Eff. fitted parameters :',handles.model_stats.nfp);
    fprintf(wfile,'%s %0d\n','Eff. degrees of freedom:',handles.model_stats.ndf);
elseif model_flag
    %If 2-distance model fit, always write <r2>/<r1> distance ratio. If nonlinear 
    %constraints exist, write their parameters. For all models, write number of
    %data observations, fitted parameters and degrees of freedom.
    fprintf(wfile,'\n%s\n','### Extra fit output ###');
    switch handles.user_model
      case {'Two_Gaussians','Two_Gaussians_hom','Two_Rice3d'}
          if handles.confit
                k1 = handles.model_pars_kcon(1);
                k2 = handles.model_pars_kcon(2);
                fprintf(wfile,'%s = %0.4f\n',handles.nlconstraint_label,handles.model_pars(k2)/handles.model_pars(k1)); 
                if handles.nconstraints > 0
                    fprintf(wfile,'%s %s\n', 'NonLin constraint :',handles.nlconstraint_label);
                    fprintf(wfile,'%s %0d\n','# NonLinCon bounds:',handles.nconstraints);
                end
          end
    end
    if isfield(handles,'model_stats') && ~isempty(handles.model_stats)
        fprintf(wfile,'%s %0d\n','Data observations :',handles.model_stats.nobs);
        fprintf(wfile,'%s %0d\n','Fitted parameters :',handles.model_stats.nfp);
        fprintf(wfile,'%s %0d\n','Degrees of freedom:',handles.model_stats.ndf);
    end
end
fclose(wfile);
handles.saved=1;
