function handback=sim_user_model(handles)
%
% Simulates distance distribution and dipolar evolution function
% for the currently selected user model to the 
% experimental dipolar evolution function
% G. Jeschke, 2006
%

full_depth=get(handles.checkbox_model_fit_depth,'Value');


set(handles.main_figure,'Pointer','watch');
drawnow;

% Set up vectors of selected fit parameters and their bounds

exflag=get(handles.exci_bandwidth_corr,'Value');

% Initialize distance axis
rmin=handles.user_fit_rmin;
rmax=handles.user_fit_rmax;
numr=round((rmax-rmin)/handles.user_fit_dr)+1;
r=linspace(rmin,rmax,numr);

sim=[];

full_par=handles.model_pars;

% Simulate fitted dipolar evolution function
if ~isempty(handles.A_tdip)
    if handles.model_mode || full_depth % for positive model mode identifiers
    my_model = handles.user_model;
    if strcmp(my_model,'Gaussian')
        my_model = 'Gaussian_and_depth';
        [sim,distr]=feval(my_model,r,handles.A_tdip,handles.model_pars);
        distr=0.01*distr/sum(distr); % normalization
    else
        distr=feval(str2func(handles.user_model),r,handles.model_pars);
        sim=get_td_fit(handles,r,distr);
    end
    else
        % Get distribution
        distr=feval(str2func(handles.user_model),r,full_par);
        distr=0.01*distr/sum(distr); % normalization

        % Simulate fitted dipolar evolution function
        if exflag
            [sim,~]=deer_sim(r,distr,handles.A_tdip,handles.bandwidth);
        else
            sim=get_td_fit(handles,r,distr);
        end;
    end;
end;

handles.model_distr=distr;
handles.model_r=r;
handles.model_sim=sim;

set(handles.main_figure,'Pointer','arrow');
drawnow;

handback=handles;
