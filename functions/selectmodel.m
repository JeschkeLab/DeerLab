%
% SELECTMODEL Optimal parametric model selection
%
%   opt = SELECTMODEL({@model1,...,@modelN},S,r,K,{'aic',...})
%   Evaluates the fits of the parametric models (model1,...,modelN) to a
%   signal (S) according to the dipolar kernel (K) and distance axis (r).
%   The models must be passed as a cell array of function handles. Each fit
%   is then evaluated according to the model selection criterions
%   ('aic','aicc','bic','rmsd') specified in the last input argument M-point
%   cell array. Function returns a M-point array containing the optimal models
%   according to each selection method.
%
%   opt = SELECTMODEL({@model1,...,@modelN},S,t,{'aic',...})
%   Evaluates the fits of the time-domain parametric models by specifying
%   the time axis (t).
%
%   opt = SELECTMODEL({@model1,...,@modelN},S,r,K,{'aic',...},{par1,...,parN})
%   The initial guess values for the parameters of each model can be passed
%   as a cell array {par1,...parN} of value vectors.
%
%   [opt,f,params] = SELECTMODEL(...)
%   Returns a cell array the method selector functionals (f) for the
%   different methods and a cell array (params) with the fitted parameters
%   for each of the evaluated models.
%
%   opt = SELECTMODEL(...,'Property',Value)
%   Additional (optional) arguments can be passed as property-value pairs.
%
%   'Lower' - Cell array containing the lower bound values for the parameters 
%             of the evaluated parametric models.
%
%   'Upper' - Cell array containing the upper bound values for the parameters 
%             of the evaluated parametric models.
%
%   See "help fitparamodel" for a detailed list of other property-value pairs
%   accepted by the function.
%

% This file is a part of DeerLab. License is MIT (see LICENSE.md).
% Copyright(c) 2019: Luis Fabregas, Stefan Stoll, Gunnar Jeschke and other contributors.


function [optima,functionals,fitparams] = selectmodel(Models,S,ax,K,Methods,param0,varargin)

if nargin<4
   error('At least four input arguments required.') 
end

isTimeDomain = false;

%Input: selectmodel(Models,S,r,K,Methods)
if ~ischar(K) && nargin<6
    param0(1:length(Models)) = {[]};
%Input: selectmodel(Models,S,t,Methods)
elseif ischar(K) && nargin<5
    param0(1:length(Models)) = {[]};
    Methods = K;
    K = [];
    isTimeDomain = true;
%Input: selectmodel(Models,S,t,Methods,'options',arg)
elseif ischar(K) && ischar(Methods)   
    varargin = [{Methods} {param0} varargin];
    param0 = {};
    param0(1:length(Models)) = {[]};
    Methods = K;
    K = [];
    isTimeDomain = true;
%Input: selectmodel(Models,S,t,Methods,param0,'options',arg)
elseif ischar(K) && ~ischar(Methods)   
    varargin = [{param0} varargin];
    param0 = Methods;
    Methods = K;
    K = [];
    isTimeDomain = true;
%Input: selectmodel(Models,S,r,K,Methods,'options',arg)
elseif ischar(param0)
    varargin = [{param0} varargin];
    param0 = {};
    param0(1:length(Models)) = {[]};
else
    if ~iscell(param0)
        error('Sixth input argument must be a cell array of inital parameter values.')
    end
end
% Input validation
%-------------------------------------------------------------------------------
if ~iscell(Models)
    error('First input must be a cell array of model functions.');
end
if ~iscell(Methods)
    Methods = {Methods};
end

%Remove the Lower and Upper options from varargin so they are not passed to fitparamodel
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

warning('off','DeerLab:parseoptional')
%Parse the optional parameters in the varargin
[Upper,Lower] = parseoptional({'Upper','Lower'},varargin2);
warning('on','DeerLab:parseoptional')

if ~isempty(Upper) && ~iscell(Upper)
    error('Upper property must be a cell array of upper bound vectors.')
end
if ~isempty(Lower) && ~iscell(Lower)
    error('Lower property must be a cell array of lower bound vectors.')
end
if length(Upper) > length(Models) || length(Lower) > length(Models)
    error('Lower/Upper bound cell array cannot exceed the number of models.')
end
%Set the bounds for the models
UpperBounds = cell(1,length(Models));
LowerBounds = cell(1,length(Models));
UpperBounds(1:length(Upper)) = Upper;
LowerBounds(1:length(Lower)) = Lower;


allowedMethodInputs = {'aic','aicc','bic','rmsd'};
for i = 1:length(Methods)
    if strcmp(Methods{i},'all')
        Methods = allowedMethodInputs;
        break;
    end
    validateattributes(Methods{i},{'char'},{'nonempty'})
    Methods{i} = validatestring(Methods{i},allowedMethodInputs);
end

% Run all parametric model fits and evaluate selection metrics
%-------------------------------------------------------------------------------
nMethods = length(Methods);
N = numel(S);
AICc = zeros(length(Models),1);
BIC = zeros(length(Models),1);
AIC = zeros(length(Models),1);
fitparams = cell(nMethods,1);
for i = 1:length(Models)
    
    if isTimeDomain
        [paramfit,Fit] = fitparamodel(S,Models{i},ax,param0{i},'Upper',UpperBounds{i},'Lower',LowerBounds{i},varargin{:});
    else
        [paramfit,Fit] = fitparamodel(S,Models{i},ax,K,param0{i},'Upper',UpperBounds{i},'Lower',LowerBounds{i},varargin{:});
    end
    fitparams{i} = paramfit;
    
    nParams = numel(paramfit);
    Q = nParams + 1;
    if isTimeDomain
        SSR = sum((S(:)-Fit(:)).^2);
    else
        SSR = sum((S(:)-K*Fit(:)).^2);
    end
    AIC(i) =  N*log(SSR/N) + 2*Q;
    AICc(i) = N*log(SSR/N) + 2*Q + 2*Q*(Q+1)/(N-Q-1);
    BIC(i) =  N*log(SSR/N) + Q*log(N);
    RMSD(i) = sqrt(1/numel(ax)*SSR);
end

% Identify optimal models based on selection criteria
%-------------------------------------------------------------------------------
optima = zeros(nMethods,1);
functionals = cell(nMethods,1);
for i = 1:nMethods
    switch Methods{i}
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

if length(functionals)==1
    functionals = functionals{1};
end
