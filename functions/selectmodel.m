%
% SELECTMODEL Optimal parametric model selection
%
%   opt = SELECTMODEL({@model1,...,@modelN},V,r,K,{'aic',...})
%   Evaluates the fits of the parametric models (model1,...,modelN) to a
%   signal (V) according to the dipolar kernel (K) and distance axis (r).
%   The models must be passed as a cell array of function handles. Each fit
%   is then evaluated according to the model selection criterions
%   ('aic','aicc','bic','rmsd') specified in the last input argument M-point
%   cell array. Function returns a M-point array containing the optimal models
%   according to each selection method.
%
%   opt = SELECTMODEL({@model1,...,@modelN},V,t,{'aic',...})
%   Evaluates the fits of the time-domain parametric models by specifying
%   the time axis (t).
%
%   opt = SELECTMODEL({@model1,...,@modelN},V,r,K,{'aic',...},{par1,...,parN})
%   The initial guess values for the parameters of each model can be passed
%   as a cell array {par1,...parN} of value vectors.
%
%   [opt,f,params,paramcis] = SELECTMODEL(...)
%   Returns a cell array the method selector functionals (f) for the
%   different methods and a cell array (params) with the fitted parameters
%   for each of the evaluated models as well as their confidence intervals (paramcis).
%
%   opt = SELECTMODEL(...,'Property',Value)
%   Additional (optional) arguments can be passed as name-value pairs.
%
%   'Lower' - Cell array containing the lower bound values for the parameters
%             of the evaluated parametric models.
%   'Upper' - Cell array containing the upper bound values for the parameters
%             of the evaluated parametric models.
%
%   See "help fitparamodel" for a detailed list of other name-value pairs
%   accepted by the function.
%

% This file is a part of DeerLab. License is MIT (see LICENSE.md).
% Copyright(c) 2019: Luis Fabregas, Stefan Stoll, Gunnar Jeschke and other contributors.


function [optima,functionals,fitparams,paramcis] = selectmodel(models,Vs,ax,Ks,methods,param0,varargin)

if nargin<4
    error('At least four input arguments required.')
end

isTimeDomain = false;

if ~ischar(Ks) && nargin<6
    % selectmodel(Models,S,r,K,Methods)
    param0(1:length(models)) = {[]};
elseif ischar(Ks) && nargin<5
    % selectmodel(Models,S,t,Methods)
    param0(1:length(models)) = {[]};
    methods = Ks;
    Ks = [];
    isTimeDomain = true;
elseif ischar(Ks) && ischar(methods)
    % selectmodel(Models,S,t,Methods,'options',arg)
    varargin = [{methods} {param0} varargin];
    param0 = {};
    param0(1:length(models)) = {[]};
    methods = Ks;
    Ks = [];
    isTimeDomain = true;
elseif ischar(Ks) && ~ischar(methods)
    % selectmodel(Models,S,t,Methods,param0,'options',arg)
    varargin = [{param0} varargin];
    param0 = methods;
    methods = Ks;
    Ks = [];
    isTimeDomain = true;
elseif ischar(param0)
    % selectmodel(Models,S,r,K,Methods,'options',arg)
    varargin = [{param0} varargin];
    param0 = {};
    param0(1:length(models)) = {[]};
else
    if ~iscell(param0)
        error('Sixth input argument must be a cell array of inital parameter values.')
    end
end

% Input validation
%-------------------------------------------------------------------------------
if ~iscell(models)
    error('First input must be a cell array of model functions.');
end
if ~iscell(methods)
    methods = {methods};
end

% Remove the Lower and Upper options from varargin so they are not passed to fitparamodel
varargin2 = [];
Idx = find(cellfun(@(x)(ischar(x) && contains(lower(x),'upper')),varargin));
varargin2 = [varargin2 varargin(Idx:Idx+1)];
varargin(Idx:Idx+1) = [];
Idx = find(cellfun(@(x)(ischar(x) && contains(lower(x),'lower')),varargin));
varargin2 = [varargin2 varargin(Idx:Idx+1)];
varargin(Idx:Idx+1) = [];
if ~isempty(varargin)
    if length(varargin)==1 && iscell(varargin{1})
        varargin = varargin{1};
    end
end

% Parse the optional parameters in the varargin
warning('off','DeerLab:parseoptional')
[Upper,Lower,GlobalWeights] = parseoptional({'Upper','Lower','GlobalWeights'},varargin2);
warning('on','DeerLab:parseoptional')

if ~isempty(Upper) && ~iscell(Upper)
    error('Upper property must be a cell array of upper bound vectors.')
end
if ~isempty(Lower) && ~iscell(Lower)
    error('Lower property must be a cell array of lower bound vectors.')
end
if length(Upper) > length(models) || length(Lower) > length(models)
    error('Lower/Upper bound cell array cannot exceed the number of models.')
end

if ~iscell(Vs)
   Vs = {Vs}; 
end
if ~isTimeDomain && ~iscell(Ks)
   Ks = {Ks}; 
end
if isTimeDomain && ~iscell(ax)
   ax = {ax}; 
end

if isempty(GlobalWeights)
    GlobalWeights = globalweights(Vs);
else
    validateattributes(GlobalWeights,{'numerical'},{'nonempty','nonnegative'},mfilename,'GlobalWeights')
    GlobalWeights = GlobalWeights/sum(GlobalWeights);
end


% Set the bounds for the models
UpperBounds = cell(1,length(models));
LowerBounds = cell(1,length(models));
UpperBounds(1:length(Upper)) = Upper;
LowerBounds(1:length(Lower)) = Lower;


allowedMethodInputs = {'aic','aicc','bic','rmsd'};
for i = 1:length(methods)
    if strcmp(methods{i},'all')
        methods = allowedMethodInputs;
        break;
    end
    validateattributes(methods{i},{'char'},{'nonempty'})
    methods{i} = validatestring(methods{i},allowedMethodInputs);
end

% Run all parametric model fits and evaluate selection metrics
%-------------------------------------------------------------------------------
nMethods = length(methods);
nModels = length(models);
AICc = zeros(nModels,1);
BIC = zeros(nModels,1);
AIC = zeros(nModels,1);
RMSD = zeros(nModels,1);
fitparams = cell(nMethods,1);
paramcis = cell(nMethods,1);
for i = 1:nModels
    if isTimeDomain
        [parfit,fit,parci] = fitparamodel(Vs,models{i},ax,param0{i},'Upper',UpperBounds{i},'Lower',LowerBounds{i},varargin{:});
    else
        [parfit,fit,parci] = fitparamodel(Vs,models{i},ax,Ks,param0{i},'Upper',UpperBounds{i},'Lower',LowerBounds{i},varargin{:});
    end
    fitparams{i} = parfit;
    paramcis{i} = parci;
    
    if isTimeDomain && ~iscell(fit)
        fit = {fit};
    end
    
    nParams = numel(parfit);
    logprob = 0;
    for idx = 1:numel(Vs)
        V = Vs{idx};
        N = numel(V);
        if isTimeDomain
            SSR = sum((V(:) - fit{idx}).^2);
        else
            K = Ks{idx};
            SSR = sum((V(:) - K*fit).^2);
        end
        Q = nParams + 1;
        logprob = logprob + GlobalWeights(idx)*N*log(SSR/N);
        RMSD(i) = RMSD(i) + GlobalWeights(idx)*sqrt(1/N*SSR);
    end
    AIC(i) =  logprob + 2*Q;
    AICc(i) = logprob + 2*Q + 2*Q*(Q+1)/(N-Q-1);
    BIC(i) =  logprob + Q*log(N);

end

% Identify optimal models based on selection criteria
%-------------------------------------------------------------------------------
optima = zeros(nMethods,1);
functionals = cell(nMethods,1);
for i = 1:nMethods
    switch methods{i}
        case 'aic'
            functional = AIC;
        case 'aicc'
            functional = AICc;
        case 'bic'
            functional = BIC;
        case 'rmsd'
            functional = RMSD;
    end
    [~,optimum] = min(functional);
    functionals{i} = functional(:);
    optima(i) = optimum;
end

if nMethods==1
    functionals = functionals{1};
end
