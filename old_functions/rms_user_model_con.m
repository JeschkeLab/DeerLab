function rms=rms_user_model_con(par,handles,r,lb,ub,exflag,mcon,kcon,ctol)
% Compute r.m.s. value of user model at given parameter vector par
%
% Modified by H.C. Hyde, 2012 
%  *Added input argument 'mcon' for distance constraints in 
%   2-distance model fits
%

% Penalize parameters out of bounds
for k=1:length(par)
    if par(k)<(lb(k)-2*ctol), rms=1.0e6; return; end;
    if par(k)>(ub(k)+2*ctol), rms=1.0e6; return; end;
end

full_par=handles.model_pars;
% Update parameters
switch handles.user_model
    case {'Two_Gaussians','Two_Rice3d'}, kmax=5;
    case 'Two_Gaussians_hom', kmax=6;
    otherwise
        error('rms_user_model_con function is only intended for minimization using fmincon')
end
for k=1:kmax
    full_par(k)=par(k);
end

t=handles.A_tdip;

if handles.model_mode % for positive model mode identifiers
    [sim,distr]=feval(str2func(handles.user_model),r,t,full_par);
    distr=0.01*distr/sum(distr); % normalization
else
    % Get distribution
    distr=feval(str2func(handles.user_model),r,full_par);
    distr=0.01*distr/sum(distr); % normalization

    % Simulate fitted dipolar evolution function
    if exflag,
        [sim,sc]=deer_sim(r,distr,handles.A_tdip,handles.bandwidth);
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
rms=sqrt(sum(diff.*diff)/(length(diff)-1));  %(HCH) fix; %OLD: rms=sqrt(sum(diff.*diff))/(length(diff)-1);

pstr=sprintf('%s%9.7f','Fitting user model. r.m.s. ',rms);
set(handles.status_line,'String',pstr);
drawnow;

% figure(13); clf;
% plot(handles.A_tdip,cluster,'k');
% hold on;
% plot(handles.A_tdip,sim,'r');
% disp(par);
% disp(rms);
% keyboard;
