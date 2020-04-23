%
% FITSIGNAL  Fit model to dipolar time-domain trace
%
%   [Vfit,Pfit,Bfit,parfit] = FITSIGNAL(V,t,r,dd,bg,ex,par0)
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
%    ex     experiment model (default @exp_4pdeer)
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
%
% Example:
%    Pfit = fitsignal(Vexp,t,r,@dd_gauss,@bg_exp,@exp_4pdeer)
%

% This file is a part of DeerLab. License is MIT (see LICENSE.md). 
% Copyright(c) 2019-2020: Luis Fabregas, Stefan Stoll and other contributors.

function [Vfit,Pfit,Bfit,parfit] = fitsignal(Vexp,t,r,dd_model,bg_model,exp_model,par0)

% Clear persistent variables in the function
clear fitsignal

if nargin<3
    error('At least three inputs (V,t,r) must be specified.');
end

validateattributes(Vexp,{'numeric'},{'vector'},mfilename,'V (1st input)');
validateattributes(t,{'numeric'},{'vector'},mfilename,'t (2nd input)');
validateattributes(r,{'numeric'},{'vector'},mfilename,'r (3rd input)');

if numel(Vexp)~=numel(t)
    error('V (1st input) and t (2nd input) must have the same number of elements.')
end

% Set defaults
if nargin<4, dd_model = 'P'; end
if nargin<5, bg_model = @bg_exp; end
if nargin<6, exp_model = @exp_4pdeer; end
if nargin<7, par0 = {[],[],[]}; end

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
if isa(exp_model,'function_handle')
    [par0_ex,lower_ex,upper_ex,N_ex] = getmodelparams(exp_model,t);
elseif ischar(exp_model) && strcmp(exp_model,'none')
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
modelidx = [ones(1,N_dd) ones(1,N_bg)*2 ones(1,N_ex)*3];
ddidx = modelidx==1;
bgidx = modelidx==2;
exidx = modelidx==3;

% Regularization settings
regtype = 'tikh';
regparam = 'aic';

if numel(par0)==0
    % Solve regularization only
    K = dipolarkernel(t,r);
    Pfit = fitregmodel(Vexp,K,r,regtype,regparam);
    Vfit = K*Pfit;
    Bfit = ones(size(Vfit));
    parfit_ = [];
else
    % Fit the parameters
    parfit_ = fitparamodel(Vexp,@Vmodel,t,par0,'Lower',lower,'Upper',upper);
  
    % Calculate the fitted signal, background, and distribution
    [Vfit,Bfit,Pfit] = Vmodel(t,parfit_);
end


% Return fitted parameter in structure
parfit.dd = parfit_(ddidx);
parfit.bg = parfit_(bgidx);
parfit.ex = parfit_(exidx);

    % General multi-pathway DEER signal model function
    function [V,B,P] = Vmodel(t,par)
        
        % Define persistents for soft-memoization of alpha-search
        persistent par_prev regparam_prev
        
        % Calculate the background and the experiment kernel matrix
        Bfcn = @(t) bg_model(t,par(bgidx));        
        if includeExperiment
            if includeBackground
                [K,B] = exp_model(t,r,par(exidx),Bfcn);
            else
                K = exp_model(t,r,par(exidx));
                B = ones(numel(t),1);
            end
        else
            K = dipolarkernel(t,r);
            if includeBackground
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
                if all(abs(par_prev-par)./par < 1e-3)
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

if contains(func2str(model),'exp_')
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
