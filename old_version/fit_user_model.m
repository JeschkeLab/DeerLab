function [handback,r,distr,sim,rms,sc] = fit_user_model(handles)
%
% Fits parameters of the currently selected user model to the 
% experimental dipolar evolution function
%
% Modified by H.C. Hyde, 2010 
%  *Added option for distance ratio constraints in 2-gaussian or 2-rice model fits
%    -implemented using DeerAnalysis parameter input window
%    -nonlinear constrained fit implemented using 'fmincon' function
%    -2 Gaussian constraint = <r2>/<r1>: min ratio = 'kmin'; max ratio = 'kmax'
%    -2 Rice constraint = <nu2>/<nu1>: min ratio = 'kmin'; max ratio = 'kmax'
%    -in general, <r2> or <nu2> is the diagonal distance, which is longer 
%     than the adjacent distance <r1> or <nu1>

maxiter = 5;
model_tolerance = 1.2; % maximum increase in rmsd compared to previous analysis that is deemed acceptable
rms0 = handles.best_rmsd; % for comparison

handles.fit_constrained = false;
confit = get(handles.checkbox_statistics,'Value');
handles.confit = confit;
if ~confit
    constrainfit = false;
end;

full_depth=get(handles.checkbox_model_fit_depth,'Value');

set(handles.main_figure,'Pointer','watch');
drawnow;

% Set up vectors of selected fit parameters and their bounds

exflag=get(handles.exci_bandwidth_corr,'Value');

handles.model_mask=zeros(1,6);
handles.model_mask(1)=get(handles.sel_par1,'Value');
handles.model_mask(2)=get(handles.sel_par2,'Value');
handles.model_mask(3)=get(handles.sel_par3,'Value');
handles.model_mask(4)=get(handles.sel_par4,'Value');
handles.model_mask(5)=get(handles.sel_par5,'Value');
handles.model_mask(6)=get(handles.sel_par6,'Value');
handles.model_mask(7)=get(handles.sel_par7,'Value');
handles.model_mask(8)=get(handles.sel_par8,'Value');

if ~confit
    switch handles.user_model
        case {'Two_Gaussians','Two_Rice3d'}
            handles.model_mask(6) = 0;
            handles.model_mask(7) = 0;
        case 'Two_Gaussians_hom'
            handles.model_mask(7) = 0;
            handles.model_mask(8) = 0;
    end;
end;
numpar=sum(handles.model_mask);
start_pars=zeros(1,numpar);
lower_bounds=zeros(1,numpar);
upper_bounds=zeros(1,numpar);
poi=0;
for k=1:8
    if handles.model_mask(k)
        poi=poi+1;
        start_pars(poi)=handles.model_pars(k);
        lower_bounds(poi)=handles.model_lower_bounds(k);
        upper_bounds(poi)=handles.model_upper_bounds(k);
    end;
end;

if confit
    % *** HC Hyde
    % If applicable, get nonlinear constraint lower/upper bounds in correct order and format
    constrainfit=0;
    mcon=struct('lb',[],'ub',[]);
    [sp,lb,ub,~]=deal([]);  %initialize start parameter values, lower bounds, upper bounds, constraint indices
    ctol=1e-6;  %fixed parameter tolerance for constrained fits
    switch handles.user_model
        case {'Two_Gaussians','Two_Rice3d'}
            np=5;  %number of TOTAL parameters for selected fit (fixed or free)
            Ip=np+[1,2];  %possible constraint parameter indices
            Ic=Ip(handles.model_mask(Ip)~=0);  %constrained parameter indices
            if ~isempty(Ic)
                if Ic==Ip(1),     mcon.lb=handles.model_pars(Ip(1));
                elseif Ic==Ip(2), mcon.ub=handles.model_pars(Ip(2));
                elseif Ic==Ip
                    mcon.lb=handles.model_pars(Ip(1));
                    mcon.ub=handles.model_pars(Ip(2));
                    if mcon.lb>mcon.ub  %check that constraints are in correct order
                        temp=sort([mcon.lb,mcon.ub]);
                        mcon.lb=temp(1);
                        mcon.ub=temp(2);
                        handles.model_pars(Ip)=deal([mcon.lb,mcon.ub]);
                    end
                end
                constrainfit=1;
                for k=1:np
                    sp(k)=handles.model_pars(k);  %start parameter
                    if handles.model_mask(k)
                        lb(k)=handles.model_lower_bounds(k);  %lower bounds
                        ub(k)=handles.model_upper_bounds(k);  %upper bounds
                    else
                        lb(k)=sp(k)-ctol;
                        ub(k)=sp(k)+ctol;
                    end
                end
            end
        case 'Two_Gaussians_hom'
            np=6;  %number of TOTAL parameters for selected fit (fixed or free)
            Ip=np+[1,2];  %possible constraint parameter indices
            Ic=Ip(handles.model_mask(Ip)~=0);  %constrained parameter indices
            if ~isempty(Ic)
                if Ic==Ip(1),     mcon.lb=handles.model_pars(Ip(1));
                elseif Ic==Ip(2), mcon.ub=handles.model_pars(Ip(2));
                elseif Ic==Ip
                    mcon.lb=handles.model_pars(Ip(1));
                    mcon.ub=handles.model_pars(Ip(2));
                    if mcon.lb>mcon.ub  %check that constraints are in correct order
                        temp=sort([mcon.lb,mcon.ub]);
                        mcon.lb=temp(1);
                        mcon.ub=temp(2);
                        handles.model_pars(Ip)=deal([mcon.lb,mcon.ub]);
                    end
                end
                constrainfit=1;
                for k=1:np
                    sp(k)=handles.model_pars(k);  %start parameter
                    if handles.model_mask(k)
                        lb(k)=handles.model_lower_bounds(k);  %lower bounds
                        ub(k)=handles.model_upper_bounds(k);  %upper bounds
                    else
                        lb(k)=sp(k)-ctol;
                        ub(k)=sp(k)+ctol;
                    end
                end
            end
    end
    handles.nconstraints = length([mcon.lb mcon.ub]);   %update number of constraints
    handles.nlconstraint_label = '';
    handles.nlconstraint_value = [];

    % Get indices of mean distance parameters if any model fit has 2 distances.
    % If nonlinear constrained fit, quit if both distances are fixed or distance parameters are invalid.
    handles.model_pars_kcon = [];
    if strcmp(handles.user_model,'Two_Gaussians') || strcmp(handles.user_model,'Two_Gaussians_hom') || strcmp(handles.user_model,'Two_Rice3d')
        dist_labels = {'<r1>','<nu1>','<mu1>';...   %list of possible mean distance labels to look for in model file
                       '<r2>','<nu2>','<mu2>'};   
        kcon=zeros(1,2);
        for j=1:size(dist_labels,1)
            [found,i]=deal(0);
            while ~found && i<size(dist_labels,2)
                i=i+1;
                temp=find(strcmpi(handles.model_pars_labels,char(dist_labels{j,i})));  %mean distance (e.g. <r1>, <r2>) parameter indices
                if ~isempty(temp)
                    kcon(j)=temp;
                    found=1;
                end
            end
        end
        if length(find(kcon))<2
            set(handles.status_line,'String','### Model file does not contain 2 valid mean distances needed for distance ratio ###');
            pause(3);
            set(handles.status_line,'String','Ready.');
            handback=handles;
            return;
        end
        k1=kcon(1);  %index of mean distance #1
        k2=kcon(2);  %index of mean distance #2
        handles.model_pars_kcon = kcon;  
        handles.nlconstraint_label = sprintf('%s/%s',char(handles.model_pars_labels{k2}),char(handles.model_pars_labels{k1}));
        if handles.model_mask(k1)==0 && handles.model_mask(k2)==0 && constrainfit==1
            switch handles.user_model
                case {'Two_Gaussians','Two_Gaussians_hom'}        
                    set(handles.status_line,'String','### Constrained 2-gaussian fit is not possible with both <r1> and <r2> fixed ###');
                case 'Two_Rice3d'
                    set(handles.status_line,'String','### Constrained 2-rice fit is not possible with both <nu1> and <nu2> fixed ###');
            end
            pause(3);
            set(handles.status_line,'String','Ready.');
            handles.fit_constrained = true;
            handback=handles;
            return;
        end
    end
    % *** end HC Hyde
end;

% Initialize distance axis
rmin=handles.user_fit_rmin;
rmax=handles.user_fit_rmax;
numr=round((rmax-rmin)/handles.user_fit_dr)+1;
r=linspace(rmin,rmax,numr);

% Actual fitting
% ************** 
if constrainfit % *** HC Hyde
    % Special Case: 2-distance model fit with constraint(s) on ratio of distribution mean distances
    opt=optimset('fmincon');
    opt=optimset(opt,'Algorithm','Active-Set');
    opt=optimset(opt,'Display','off',...
                     'MaxIter',300,...
                     'MaxFunEvals',3000,...
                     'TolFun',1e-10,...
                     'TolCon',1e-10,...
                     'DiffMinChange',1e-8,...
                     'DiffMaxChange',0.1);
    % Adjust initial time constant values relative to constraints (if necessary)
    % If either mcon.lb or mcon.ub is empty, the corresponding 'if' statement will be ignored.
    kcurrent = sp(k2)/sp(k1);
    if kcurrent<mcon.lb || kcurrent>mcon.ub
        if handles.model_mask(k1)==1 && handles.model_mask(k2)==1  %<r1> fitted, <r2> fitted
            if kcurrent<mcon.lb, ss=-1; end
            if kcurrent>mcon.ub, ss=+1; end
            if ~isempty(mcon.ub)
                 aa=abs(kcurrent-mcon.ub)/kcurrent;
            else
                aa=abs(kcurrent-mcon.lb)/kcurrent;
            end
            if ~isempty(aa)
                sp(k1)=(1+ss*aa/2)*sp(k1);
                sp(k2)=(1-ss*aa/2)*sp(k2);
            end 
        elseif handles.model_mask(k1)==1 && handles.model_mask(k2)==0  %<r1> fitted, <r2> fixed
            if kcurrent<mcon.lb, sp(k1)=sp(k2)/mcon.lb; end
            if kcurrent>mcon.ub, sp(k1)=sp(k2)/mcon.ub; end
        elseif handles.model_mask(k1)==0 && handles.model_mask(k2)==1  %<r1> fixed, <r2> fitted
            if kcurrent<mcon.lb, sp(k2)=sp(k1)*mcon.lb; end
            if kcurrent>mcon.ub, sp(k2)=sp(k1)*mcon.ub; end
        end
    end
    rms = 2*rms0;
    iter = 0;
    while rms > model_tolerance*rms0 && iter < maxiter
        iter = iter + 1;
        [fit_pars,rms] = fmincon(@rms_user_model_con,sp,[],[],[],[],lb,ub,@mean_ratio_con,opt,handles,r,lb,ub,exflag,mcon,kcon,ctol);
        start_pars = fit_pars;
    end
else  
    % General Case: All other models with general lower/upper bounds
    rms = 2*rms0;
    iter = 0;
    while rms > model_tolerance*rms0 && iter < maxiter
        iter = iter + 1;
        [fit_pars,rms] = fminsearch(@rms_user_model,start_pars,[],handles,r,lower_bounds,upper_bounds,exflag);
        start_pars = fit_pars;
    end
end

if ~handles.user_model_trial
    if rms > model_tolerance*rms0
        set(handles.status_line,'String','### Parametrized model fits significantly worse than best previous analysis ###');
    else
        set(handles.status_line,'String','Ready.');
    end
end
% Update parameters
% Check for constrained fit selection; in such case update all model parameters (but not constraints).
% Otherwise, only update fitted parameters.
if confit && handles.nconstraints>0 % *** HC Hyde
    switch handles.user_model
      case {'Two_Gaussians','Two_Rice3d'}
        kmax=5;
        for k=1:kmax, handles.model_pars(k)=fit_pars(k); end
      case 'Two_Gaussians_hom'
        kmax=6;
        for k=1:kmax, handles.model_pars(k)=fit_pars(k); end
    end
else
    poi=0;
    kmax=8;
    for k=1:kmax
        if handles.model_mask(k)
            poi=poi+1;
            handles.model_pars(k)=fit_pars(poi);
        end
    end    
end

% Update displayed values
pstr=sprintf('%0.4g',handles.model_pars(1));
set(handles.par1_edit,'String',pstr);
set(handles.par1_edit,'ForegroundColor','b');
pstr=sprintf('%0.4g',handles.model_pars(2));
set(handles.par2_edit,'String',pstr);
set(handles.par2_edit,'ForegroundColor','b');
pstr=sprintf('%0.4g',handles.model_pars(3));
set(handles.par3_edit,'String',pstr);
set(handles.par3_edit,'ForegroundColor','b');
pstr=sprintf('%0.4g',handles.model_pars(4));
set(handles.par4_edit,'String',pstr);
set(handles.par4_edit,'ForegroundColor','b');
pstr=sprintf('%0.4g',handles.model_pars(5));
set(handles.par5_edit,'String',pstr);
set(handles.par5_edit,'ForegroundColor','b');
pstr=sprintf('%0.4g',handles.model_pars(6));
set(handles.par6_edit,'String',pstr);
set(handles.par6_edit,'ForegroundColor','b');
pstr=sprintf('%0.4g',handles.model_pars(7));
set(handles.par7_edit,'String',pstr);
set(handles.par7_edit,'ForegroundColor','b');
pstr=sprintf('%0.4g',handles.model_pars(8));
set(handles.par8_edit,'String',pstr);
set(handles.par8_edit,'ForegroundColor','b');


% Simulate fitted dipolar evolution function
sc=1;
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
    % Get distance distribution for fitted parameters
    distr=feval(str2func(handles.user_model),r,handles.model_pars);
    if exflag
        [sim,sc]=deer_sim(r,distr,handles.A_tdip,handles.bandwidth);
    else
        sim=get_td_fit(handles,r,distr);
    end;
end;

if ~handles.user_model_trial
    handles.moddepth_suppression=sc;
    handles.A_sim=sim;
    handles.A_r=r;
    handles.mask=ones(size(distr));
    handles.A_distr=distr;
    handles.A_low=distr;
    handles.A_high=distr;
    handles.updated=1;
    pstr=sprintf('%8.6f',rms);
    set(handles.distr_rms,'String',pstr);
end

set(handles.main_figure,'Pointer','arrow');
drawnow;

handback=handles;
