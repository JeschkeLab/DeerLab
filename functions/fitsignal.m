%
% FITSIGNAL  Fit model to dipolar time-domain trace
%
%   [Pfit,Vfit,Bfit,parfit] = FITSIGNAL(V,t,r,dd,bg,ex,par0)
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
%    dd     function handle to distribution model (for parametric distribution)
%           or [] for parameter-free distribution; default: []
%    bg     function handle to background model, or [] if no background should
%           be included; default []
%    ex     function handle to experiment model; default []
%    par0   starting parameters, 3-element cell array {par0_dd,par0_bd,par0_ex}
%           default: {[],[],[]} (automatic choice)
%
%  Output:
%    Pfit   fitted distance distribution
%    Vfit   fitted time-domain signal
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

function [Pfit,Vfit,Bfit,parfit] = fitsignal(Vexp,t,r,dd_model,bg_model,ex_model,par0)

if nargin<3
    error('At least three inputs (V,t,r) must be specified.');
end

validateattributes(Vexp,{'numeric'},{'vector'},mfilename,'V (1st input)');
validateattributes(t,{'numeric'},{'vector'},mfilename,'t (2nd input)');
validateattributes(r,{'numeric'},{'vector'},mfilename,'r (3rd input)');

if numel(Vexp)~=numel(t)
    error('V (1st input) and t (2nd input) must have the same number of elements.')
end

if nargin<4, dd_model = []; end
if nargin<5, bg_model = []; end
if nargin<6, ex_model = []; end
if nargin<7, par0 = {[],[],[]}; end

if ~isempty(par0)
    if ~iscell(par0) || numel(par0)~=3
        error('Initial parameters (7th input) must be a 3-element cell array.')
    end
end

% Get information about distance distribution parameters
if isa(dd_model,'function_handle')
    [par0_dd,lower_dd,upper_dd,N_dd] = getmodelparams(dd_model);
elseif isempty(dd_model)
    par0_dd = [];
    lower_dd = [];
    upper_dd = [];
    N_dd = 0;
else
    error('Distribution model (4th input) must either be a function handle (for a parametric model) or [] (for a parameter-free distribution).')
end

% Get information about background parameters
if isa(bg_model,'function_handle')
    [par0_bg,lower_bg,upper_bg,N_bg] = getmodelparams(bg_model);
elseif isempty(bg_model)
    par0_bg = [];
    lower_bg = [];
    upper_bg = [];
    N_bg = 0;
else
    error('Background model (5th input) must either be a function handle, or [] if no background should be fitted.')
end

% Get information about experiment parameters
if isa(ex_model,'function_handle')
    [par0_ex,lower_ex,upper_ex,N_ex] = getmodelparams(ex_model);
elseif isempty(ex_model)
    par0_ex = [];
    lower_ex = [];
    upper_ex = [];
    N_ex = 0;
else
    error('Experiment model (6th input) must either be a function handle, or [] if no experimental parameters should be fitted.')
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
        
        % Calculate the background and the experiment kernel matrix
        if isa(bg_model,'function_handle')
            Bfcn = @(t) bg_model(t,par(bgidx));
        else
            Bfcn = [];
        end
        if isa(ex_model,'function_handle')
            [K,B] = ex_model(t,r,par(exidx),Bfcn);
        else
            K = dipolarkernel(t,r);
            if ~isempty(Bfcn)
                B = Bfcn(t);
            else
                B = ones(size(t));
            end
        end
        
        % Get the distance distribution from the model or via regularization
        if isa(dd_model,'function_handle')
            P = dd_model(r,par(ddidx));
        else
            P = fitregmodel(Vexp,K,r,regtype,regparam);
        end
        
        % Compute the dipolar signal
        V = K*P;
        
    end
    
end

function [par0,lo,up,N] = getmodelparams(model)

info = model();
par0 = [info.parameters.default];
range = [info.parameters.range];
lo = range(1:2:end-1);
up = range(2:2:end);
N = numel(par0);

end
