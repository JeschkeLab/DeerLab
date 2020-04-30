%
% FITSIGNAL  Fit model to dipolar time-domain trace
%
%   [Vfit,Pfit,Bfit,parfit,parci,stats] = FITSIGNAL(V,t,r,dd,bg,ex,par0)
%   __ = FITSIGNAL(V,t,r,dd,bg)
%   __ = FITSIGNAL(V,t,r,dd)
%   __ = FITSIGNAL(V,t,r)
%
%   Fits a dipolar model to the experimental signal in V, using distance axis r.
%   The model is specified by the distance distribution (dd), the background
%   (bg), and the experiment (ex).
%
%   FITSIGNAL can handle both parametric and non-parametric distance
%   distribution models.
%
%  Input:
%    V      time-domain signal to fit (N-element vector)
%    t      time axis, in microseconds (N-element vector)
%    r      distance axis, in nanometers (M-element vector)
%    dd     distance distribution model (default 'P')
%           - function handle to parametric distribution model
%           - 'P' to indicate parameter-free distribution (default)
%           - 'none' to indicate no distribution, i.e. only background
%    bg     background model (default @bg_exp)
%           - function handle to parametric background model
%           - 'none' to indicate no background decay
%    ex     experiment model (default @ex_4pdeer)
%           - function handle to experiment model
%           - 'none' to indicate simple dipolar oscillation (mod.depth = 1)
%    par0   starting parameters, 3-element cell array {par0_dd,par0_bd,par0_ex}
%           default: {[],[],[]} (automatic choice)
%
%  Output:
%    Vfit   fitted time-domain signal
%    Pfit   fitted distance distribution
%    Bfit   fitted background decay
%    parfit structure with fitted parameters
%           .dd  fitted parameters for distance distribution model
%           .bg  fitted parameters for background model
%           .ex  fitted parameters for experiment model
%    parci structure with confidence intervals for parameter, similar to parfit
%    stats goodness of fit statistical estimators, N-element structure array 

% Example:
%    Vfit = fitsignal(Vexp,t,r,@dd_gauss,@bg_exp,@ex_4pdeer)
%

% This file is a part of DeerLab. License is MIT (see LICENSE.md). 
% Copyright(c) 2019-2020: Luis Fabregas, Stefan Stoll and other contributors.

function [Vfit,Pfit,Bfit,parfit,parci,stats] = fitsignal(Vexp,t,r,dd_model,bg_model,ex_model,par0)

if nargin<3
    error('At least three inputs (V,t,r) must be specified.');
end

validateattributes(Vexp,{'numeric'},{'vector'},mfilename,'V (1st input)');
validateattributes(t,{'numeric'},{'vector'},mfilename,'t (2nd input)');
validateattributes(r,{'numeric'},{'vector'},mfilename,'r (3rd input)');

if numel(Vexp)~=numel(t)
    error('V (1st input) and t (2nd input) must have the same number of elements.')
end

% Regularization settings
regtype = 'tikh';
regparam = 'aic';
alphaOptThreshold = 1e-3; % relative parameter change threshold for reoptimizing alpha


% Set defaults
if nargin<4, dd_model = 'P'; end
if nargin<5, bg_model = @bg_hom3d; end
if nargin<6, ex_model = @ex_4pdeer; end
if nargin<7, par0 = {[],[],[]}; end

calculateCI = nargout>=5 || nargout==0;

if ~isempty(par0)
    if ~iscell(par0) || numel(par0)~=3
        error('Initial parameters (7th input) must be a 3-element cell array.')
    end
end

% Get information about distance distribution parameters
par0_dd = [];
lower_dd = [];
upper_dd = [];
N_dd = 0;
includeForeground = true;
if isa(dd_model,'function_handle')
    [par0_dd,lower_dd,upper_dd,N_dd] = getmodelparams(dd_model);
    parfreeDistribution = false;
elseif ischar(dd_model) && strcmp(dd_model,'P')
    parfreeDistribution = true;
elseif ischar(dd_model) && strcmp(dd_model,'none')
    includeForeground = false;
    parfreeDistribution = false;
else
    error('Distribution model (4th input) must either be a function handle, ''P'', or ''none''.')
end

% Get information about background parameters
par0_bg = [];
lower_bg = [];
upper_bg = [];
N_bg = 0;
includeBackground = true;
if isa(bg_model,'function_handle')
    [par0_bg,lower_bg,upper_bg,N_bg] = getmodelparams(bg_model);
elseif ischar(bg_model) && strcmp(bg_model,'none')
    includeBackground = false;
else
    error('Background model (5th input) must either be a function handle, or ''none''.')
end

% Get information about experiment parameters
par0_ex = [];
lower_ex = [];
upper_ex = [];
N_ex = 0;
includeExperiment = true;
if isa(ex_model,'function_handle')
    [par0_ex,lower_ex,upper_ex,N_ex] = getmodelparams(ex_model,t);
elseif ischar(ex_model) && strcmp(ex_model,'none')
    includeExperiment = false;
else
    error('Experiment model (6th input) must either be a function handle, or ''none''.')
end

% Catch nonsensical situation
if ~includeForeground && ~includeBackground
    error('Cannot fit anything without distribution model and without background model.')
end

% Combine all parameters into a single vector
if isempty(par0{1}), par0{1} = par0_dd; end
if isempty(par0{2}), par0{2} = par0_bg; end
if isempty(par0{3}), par0{3} = par0_ex; end
par0 = [par0{1} par0{2} par0{3}];
 
lower = [lower_dd lower_bg lower_ex];
upper = [upper_dd upper_bg upper_ex];

% Build index vectors for accessing parameter subsets
modelidx = [ones(1,N_dd) ones(1,N_bg)*2 ones(1,N_ex)*3].';
ddidx = modelidx==1;
bgidx = modelidx==2;
exidx = modelidx==3;

if numel(par0)==0
    % Solve regularization only
    K = dipolarkernel(t,r);
    Pfit = fitregmodel(Vexp,K,r,regtype,regparam);
    Vfit = K*Pfit;
    Bfit = ones(size(Vfit));
    parfit_ = [];
    parci_ = [];
else
    % Keep track of alpha and parameter vector across iterations, to avoid
    % doing alpha optimizations if parameter vector doesn't change much
    par_prev = [];
    regparam_prev = [];
    % Fit the parameters
    if calculateCI
        [parfit_,~,parci_,~,stats] = fitparamodel(Vexp,@Vmodel,t,par0,'Lower',lower,'Upper',upper,'TolFun',1e-5);
    else
        parfit_ = fitparamodel(Vexp,@Vmodel,t,par0,'Lower',lower,'Upper',upper,'TolFun',1e-5);
    end
  
    % Calculate the fitted signal, background, and distribution
    alpha = regparam; % use original setting for final run
    [Vfit,Bfit,Pfit] = Vmodel(t,parfit_);
end

% Return fitted parameter in structure
parfit_ = parfit_(:);
parfit.dd = parfit_(ddidx);
parfit.bg = parfit_(bgidx);
parfit.ex = parfit_(exidx);

if calculateCI
    parci.dd = parci_(ddidx,:);
    parci.bg = parci_(bgidx,:);
    parci.ex = parci_(exidx,:);
end

% Plotting
if nargout==0
    subplot(2,1,1);
    plot(t,Vexp,t,Vfit)
    axis tight
    grid on
    xlabel('time (us)');
    ylabel('V');
    legend('exp','fit');
    subplot(2,1,2);
    plot(r,Pfit);
    xlabel('distance (nm)');
    axis tight
    ylabel('P (nm^{-1})');
    grid on
    
    disp('Fitted parameters and confidence intervals')
    str = '  %s(%d):   %10f  (%10f, %10f)  %s (%s)\n';
    if numel(parfit.dd)>0
        pars = dd_model().parameters;
        for p = 1:numel(parfit.dd)
            fprintf(str,'dd',p,parfit.dd(p),...
                parci.dd(p,1),parci.dd(p,2),pars(p).name,pars(p).units);
        end
    end
    if numel(parfit.bg)>0
        pars = bg_model().parameters;
        for p = 1:numel(parfit.bg)
            fprintf(str,'bg',p,parfit.bg(p),...
                parci.bg(p,1),parci.bg(p,2),pars(p).name,pars(p).units)
        end
    end
    if numel(parfit.ex)>0
        pars = ex_model(t).parameters;
        for p = 1:numel(parfit.ex)
            fprintf(str,'ex',p,parfit.ex(p),...
                parci.ex(p,1),parci.ex(p,2),pars(p).name,pars(p).units)
        end
    end
end

    % General multi-pathway DEER signal model function
    function [V,B,P] = Vmodel(t,par)
        % Calculate the background and the experiment kernel matrix
        if includeExperiment
            pathinfo = ex_model(t,par(exidx));
            if includeBackground
                Bfcn = @(t,lam) bg_model(t,par(bgidx),lam);
                B = dipolarbackground(t,pathinfo,Bfcn);
            else
                Bfcn = [];
                B = ones(numel(t),1);
            end
            K = dipolarkernel(t,r,pathinfo,Bfcn);
        else
            K = dipolarkernel(t,r);
            if includeBackground
                Bfcn = @(t) bg_model(t,par(bgidx));
                B = Bfcn(t);
            else
                B = ones(numel(t),1);
            end
        end
        
        % Get the distance distribution
        if includeForeground
            
            % Use the alpha-search settings by default
            alpha = regparam;
            % If the parameter vectors has not changed by much...
            if ~isempty(par_prev)
                if all(abs(par_prev-par)./par < alphaOptThreshold)
                    % ...use the alpha optimized in the previous iteration
                    alpha = regparam_prev;
                end
            end
            par_prev = par;
                        
            if parfreeDistribution
                [P,~,regparam_prev] = fitregmodel(Vexp,K,r,regtype,alpha);
            else
                P = dd_model(r,par(ddidx));
            end
        else
            P = zeros(numel(t),1);
        end
        
        % Compute the total signal
        if includeForeground
            V = K*P;
        else
            V = B;
        end
    end
    
end

function [par0,lo,up,N] = getmodelparams(model,t)

if contains(func2str(model),'ex_')
    info = model(t);
else
    info = model();
end
par0 = [info.parameters.default];
range = [info.parameters.range];
lo = range(1:2:end-1);
up = range(2:2:end);
N = numel(par0);

end
