function save_list_result(handles,list_name,r,distr,tstd,mean_dipevo,meanr,stdr,mean_sig,stdstd,mean_m3,stdm3,arr_mod_norm,cmpmat),

%
% Saves the distance distribution, the experimental and simulated dipolar
% evolution functions, and the result of the moment analysis
%
% method    string that is used as middle part of the file name
%
% Content and file name convention:
% - results are saved in the directory specified in handles.project
% - filenames are derived from the name of the input dataset
%   (handles.bas_name)
% - the method string is appended to the name of the data set
% - five ASCII output files are generated, which are distinguished by another
%   appendix to the file name:
%   _distr  distance distribution, two columns, distance axis and
%           distribution
%   _res    main parameters of the analysis and moment analysis of the
%           distribution
%   _mean   mean dipolar evolution function (S/N-weighted)
%   _cmp    modulation depth normalized primary data
%   _diff   matrix of differences of modulation depth normalized data
% - files _distr, _cmp, and _diff have the extension '.dat', file _res has
%   the extension '.txt'
% - time axes are in microseconds, distance axes are in nanometers
%

fname_distr=[list_name '_distr.dat'];
fname_res=[list_name '_res.txt'];
fname_cmp=[list_name '_cmp.dat'];
fname_mean=[list_name '_mean.dat'];
fname_diff=[list_name '_diff.dat'];

tdip=tstd;
dipevo=mean_dipevo;

data1=[r' distr'];
save(fname_distr,'data1','-ascii');
data3=[tdip' dipevo'];
save(fname_mean,'data3','-ascii');
save(fname_cmp,'arr_mod_norm','-ascii');
save(fname_diff,'cmpmat','-ascii');

% Determine analysis method
method='Approximate Pake transformation';
APT_flag=get(handles.select_APT,'Value');
Tikh_flag=get(handles.select_Tikhonov,'Value');
if Tikh_flag,
    method='Tikhonov regularization';
	Lcurve_flag=get(handles.select_L_curve,'Value');
    if Lcurve_flag,
        method=[method ' with L curve computation'];
    end;
end;
model_flag=get(handles.select_model,'Value');
if model_flag,
    method=['Fit by model: ' handles.user_model];
end;

% Determine background model
bckg_model='No';
noflag=1;
pflag=get(handles.bckg_poly,'Value');
if pflag, bckg_model='Polynomial'; end;
hflag=get(handles.bckg_homogeneous,'Value');
if hflag,
    bckg_model='Homogeneous ';
    pstr=num2str(handles.hom_dim);
    bckg_model=[bckg_model pstr ' dimensional'];
    dflag=get(handles.bckg_fit_dim,'Value');
end;
uflag=get(handles.bckg_exp,'Value');
if uflag, 
    bckg_model='Fixed polynomial (experimental)';
end;

wfile=fopen(fname_res,'w+');
fprintf(wfile,'%s%s%s\n\n','>>> DEER analysis of list of data sets ',list_name,' <<<');
fprintf(wfile,'%s%s%s\n\n','- ',method,' -');
fprintf(wfile,'\n%s\n','### Pre-processing ###');
fprintf(wfile,'%s%d%s\n','Time shift: ',handles.zerotime,' ns');
fprintf(wfile,'%s%d%s\n','Cutoff at : ',handles.cutoff,' ns');
fprintf(wfile,'%s%s\n',bckg_model,' background correction');
fprintf(wfile,'%s%d%s\n','Starting at : ',handles.bckg_start,' ns');
if pflag | uflag,
	fprintf(wfile,'%s\n','Background polynomial (for logarithmic data):');
    fprintf(wfile,'%s','(');
    for k=1:length(handles.polynomial),
        fprintf(wfile,'%0.4g',handles.polynomial(k));
        if k<length(handles.polynomial),
            fprintf(wfile,'%s',', ');
        else,
            fprintf(wfile,'%s\n',')');
        end;
    end;
end;
lp_flag=get(handles.long_pass,'Value');
if lp_flag,
    fprintf(wfile,'%s%5.2f%s\n','Distance long-pass filter at ',handles.longpass_min,' nm');
end;
if APT_flag,
    fprintf(wfile,'%s%5.2f%s\n','Distance-domain smoothing filter width: ',handles.DDS,' nm');
end;
if Tikh_flag | model_flag,
    exflag=get(handles.exci_bandwidth_corr,'Value');
    if exflag,
        fprintf(wfile,'%s%5.2f%s\n','Correction for excitation bandwidth of ',handles.bandwidth,' MHz');
    end;
end;
fprintf(wfile,'\n%s\n','### Global distance distribution ###');
fprintf(wfile,'%s%5.2f%s%5.2f%s\n','Mean distance  <r> ',meanr,' nm, with std. dev. ',stdr,' nm');
fprintf(wfile,'%s%5.2f%s%5.2f%s\n','Std. dev. sigma(r) ',mean_sig,' nm, with std. dev. ',stdstd,' nm');
fprintf(wfile,'%s%5.2f%s%5.2f%s\n','3rd moment         ',mean_m3,' nm, with std. dev. ',stdm3,' nm');
fclose(wfile);