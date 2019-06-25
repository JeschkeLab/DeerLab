function rms=rms_user_model(par,handles,r,lower_bounds,upper_bounds,exflag)
% Compute r.m.s. value of user model at given parameter vector par
%

% Penalize parameters out of bounds
for k=1:length(par)
    if par(k)<lower_bounds(k), rms=1.0e6; return; end;
    if par(k)>upper_bounds(k), rms=1.0e6; return; end;
end;

full_par=handles.model_pars;
% Update parameters
poi=0;
for k=1:8
    if handles.model_mask(k)
        poi=poi+1;
        full_par(k)=par(poi);
    end;
end;

t=handles.A_tdip;

if handles.model_mode % for positive model mode identifiers
    [sim,distr]=feval(str2func(handles.user_model),r,t,full_par);
    distr=0.01*distr/sum(distr); % normalization
else
    % Get distribution
    distr=feval(str2func(handles.user_model),r,full_par);
    distr=0.01*distr/sum(distr); % normalization

    % Simulate fitted dipolar evolution function
    if exflag
        sim = deer_sim(r,distr,handles.A_tdip,handles.bandwidth);
    else
        sim=get_td_fit(handles,r,distr);
    end;
end;


cluster=handles.A_cluster;

modsim=ones(size(sim))-sim;
modexp=ones(size(cluster))-cluster;
sc=sum(modexp.*modexp)/sum(modsim.*modexp);
if sc<0, rms=1e6; return; end;
sim=ones(size(modsim))-sc*modsim;
diff=sim-cluster;
rms=sqrt(sum(diff.*diff)/(length(diff)-1));  %(HCH) fix;  %OLD: rms=sqrt(sum(diff.*diff))/(length(diff)-1);

if ~handles.user_model_trial
    pstr=sprintf('%s%9.6f','Fitting user model. r.m.s. ',rms);
    set(handles.status_line,'String',pstr);
    drawnow;
end

% figure(13); clf;
% plot(handles.A_tdip,cluster,'k');
% hold on;
% plot(handles.A_tdip,sim,'r');
% disp(par);
% disp(rms);
% keyboard;
