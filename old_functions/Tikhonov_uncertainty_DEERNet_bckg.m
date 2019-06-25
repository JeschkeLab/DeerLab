function handles = Tikhonov_uncertainty_DEERNet_bckg(handles)
% Computes a Tikhonov regularization uncertainty estimate based on the
% ensemble of background models from DEERNet, noise level in the Monte
% Carlo procedure is twice the estimated noise level (adding as much noise
% as is already there)
%
% 5 pseudo-random noise instances per background model are tested
% background models are considered to be independent of the additional noise
%
% Tikhonov_uncertainty_DEERNet_bckg(handles)
%
% G. Jeschke, 31.3.2018

num_noise = 5;
nfac = 0.5; % noise level 

tdip0 = handles.A_tdip;
dipevo0 = handles.A_dipevo;
cluster0 = handles.A_cluster;

ghost_suppression=get(handles.checkbox_ghost,'Value');

bckg_ensemble = handles.deernet_ensemble_bckg;
[mb,~] = size(bckg_ensemble);

nt = mb*num_noise;

trial_dipevo = zeros(nt,length(dipevo0));
trial_cluster = zeros(nt,length(cluster0));
trial_distr = zeros(nt,4196);
trial_sim = zeros(nt,length(tdip0));
rms_vec= zeros(1,nt);
mod_depths = zeros(1,nt);

set(handles.status_line,'String','Computing validation with DEERNet background ensemble...');

tic;
poi = 0;
for kb = 1:mb
    for kn = 1:num_noise
        noise = nfac*handles.best_rmsd*randn(size(handles.A_deernet_vexp));
        cluster = (handles.A_deernet_vexp+noise)./bckg_ensemble(kb,:);
        cluster = cluster/max(cluster);
        my_offset = 1-handles.A_depth;
        if ghost_suppression
            cluster=cluster.^(1/(handles.spins_per_object-1));
            my_offset = my_offset^(1/(handles.spins_per_object-1)); 
        end
        dipevo = cluster - my_offset;
        handles.A_dipevo = dipevo/max(dipevo);
        handles.A_tdip = handles.A_deernet_t;
        handles.A_cluster = cluster;
        axes(handles.original_data);
        cla;
        set(gca,'FontSize',8);
        plot(handles.A_deernet_t,handles.A_deernet_vexp+noise,'k');
        hold on;
        plot(handles.A_deernet_t,bckg_ensemble(kb,:),'r');
        axes(handles.dipolar_evolution);
        cla;
        set(gca,'FontSize',8);
        plot(handles.A_deernet_t,cluster,'k');
        [A_r,distr,sim,rms,sc]=get_Tikhonov(handles);
        axes(handles.distance_distribution);
        set(handles.distance_distribution,'NextPlot','replace');
        cla;
        plot(A_r,distr,'k');
        set(gca,'FontSize',8);
        set(handles.status_line,'String',sprintf('Trial %i/%i finished.',(kb-1)*num_noise+kn,nt));
        drawnow
        if ~isempty(distr)
            poi=poi+1;
            mod_depths(poi) = 1 - bckg_ensemble(kb,1);
            modsim=ones(size(sim))-sim;
            sim1=ones(size(modsim))-sc*modsim;
            trial_dipevo(poi,:)=dipevo;
            trial_cluster(poi,:)=cluster;
            trial_distr(poi,1:length(distr))=distr;
            trial_sim(poi,:)=sim1;
            rms_vec(poi)=rms;
        end;
    end
end
trial_dipevo = trial_dipevo(1:poi,:);
trial_cluster = trial_cluster(1:poi,:);
trial_distr = trial_distr(1:poi,1:length(distr));
trial_sim = trial_sim(1:poi,:);
rms_vec = rms_vec(1:poi);
handles.successful_trials=poi;
endtime=toc;
set(handles.status_line,'String',sprintf('Validation finished after %8.0f s',endtime));
distr_std=std(trial_distr);
mean_distr=mean(trial_distr);
mean_distr(mean_distr<0) = 0;
handles.mean_distr=mean_distr;
handles.A_distr = mean_distr;
handles.distr_std=distr_std;
mean_depth=mean(mod_depths);
std_depth=std(mod_depths);
handles.mean_depth=mean_depth;
handles.std_depth=std_depth;

handles.trial_distr=trial_distr;
handles.trial_sim=trial_sim;
handles.trial_cluster=trial_cluster;
handles.trial_dipevo=trial_dipevo;
handles.trial_rmsd=rms_vec;

dlow = min(trial_distr);
dlow(dlow<0) = 0;
dhigh = max(trial_distr);
handles.A_low=dlow;
handles.A_high=dhigh;

handles.A_tdip=tdip0;
handles.A_dipevo=dipevo0;
handles.A_cluster=cluster0;

if poi > 1
    set(handles.error_estimate,'Value',1);
end

function [r,distr,sim,rms,sc]=get_Tikhonov(handles)
%
% Tikhonov regularization, with or without excitation bandwidth correction
% and with or without computation of a whole L curve
% Matlab-internal version
%
% G. Jeschke, 2015-2018
%


[r,distr] = get_Tikhonov_new(handles,handles.regpar);
exflag=get(handles.exci_bandwidth_corr,'Value');
if exflag
    sim = deer_sim(r,distr,handles.A_tdip,handles.bandwidth);
else
    sim = get_td_fit(handles,r,distr);
end;
sim = real(sim);
modsim=ones(size(sim))-sim;
modexp=ones(size(handles.A_cluster))-handles.A_cluster;
sc=sum(modexp.*modexp)/sum(modsim.*modexp);
sc = real(sc);
sim1=ones(size(modsim))-sc*modsim;
diff0=handles.A_cluster-sim1;
rms=sqrt(sum(diff0.*diff0)/(length(diff0)-1));

