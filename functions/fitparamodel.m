%
% FITPARAMODEL Fits a time- or distance-domain parametric model to one (or several) signals
%
%   [param,Vfit,paramci,modelci,stats] = FITPARAMODEL(V,@model,t)
%   [param,Vfit,paramci,modelci,stats] = FITPARAMODEL(V,@model,r,K)
%   Fits the N-point signal (V) to a M-point parametric model (@model) given an
%   M-point distance/time axis (r/t). For distance-domain fitting, provide
%   the NxM point kernel matrix (K). The fitted model corresponds to a parametric model
%   calculated by the passed function handle (@model). The fitted parameters (param)
%   are returned as the first output argument, their 95% confidence intervals (paramci) are
%   returned as the third output, the fitted model as the second output, and the
%   corresponding 95% confidence bands (modelci) as the fourth output. If
%   more than one confidence level is requested, (paramci) and (modelci)
%   are given as cell arrays containing the confidence intervals at the different
%   confidence levels. A structure containing different statistical
%   estimators of goodness of fit is returned as (stats).
%
%   [param,Vfit,paramci,modelci] = FITPARAMODEL(V,@model,t,param0)
%   [param,Vfit,paramci,modelci] = FITPARAMODEL(V,@model,r,K,param0)
%   The initial guess of the model parameters can be passed as a last
%   argument (param0). This is optional for DeerLab model functions. If (@model)
%   is a user-defined function handle, (param0) is required.
%
%   [param,Vfit,paramci,modelci] = FITPARAMODEL({V1,V2,___},@model,{t1,t2,___},param0)
%   [param,Vfit,paramci,modelci] = FITPARAMODEL({V1,V2,___},@model,r,{K1,K2,___},param0)
%   Pass multiple signals/kernels to enable global fitting of a single parametric
%   model to all data. The global fit weights are automatically computed according
%   to their contribution to ill-posedness.
%
%   [param,Vfit,paramci,modelci] = FITPARAMODEL(___,'Property',Values)
%   Additional (optional) arguments can be passed as name-value pairs.
%
% The properties to be passed as options can be set in any order.
%
%   'Solver' - Solver to be used to solve the minimization problems
%           'lsqnonlin' - Non-linear constrained least-squares (toolbox)
%             'fmincon' - Non-linear constrained minimization (toolbox)
%             'nlsqbnd' - Non-linear constrained least-squares (free)
%       'fminsearchcon' - Non-linear constrained minimization (free)
%          'fminsearch' - Unconstrained minimization
%   'GlobalWeights' - Array of weighting coefficients for the individual signals in
%                     global fitting.
%   'Algorithm' - Algorithm to be used by the solvers (see fmincon or
%                 lsqnonlin documentation)
%   'Lower' - Lower bound on the model parameters
%   'Upper' - Upper bound on the model parameters
%   'TolFun' - Optimizer function tolerance
%   'MaxIter' - Maximum number of optimizer iterations
%   'MaxFunEvals' - Maximum number of optimizer function evaluations
%   'MultiStart' - Number of starting points for global optimization
%   'ConfidenceLevel' - Confidence level(s) for confidence intervals
%   'Verbose' - Display options for the solvers:
%                 'off' - no information displayed
%                 'final' - display solver exit message
%                 'iter-detailed' - display state of solver at each iteration
%               See MATLAB doc optimoptions for detailed explanation
%

% This file is a part of DeerLab. License is MIT (see LICENSE.md).
% Copyright(c) 2019-2020: Luis Fabregas, Stefan Stoll and other contributors.

function [parfit,modelfit,parci,modelci,stats] = fitparamodel(V,model,ax,K,StartParameters,varargin)

% Input parsing & validation
%-------------------------------------------------------------------------------
if nargin<3
    error('At least three inputs (V, model, t) are required.');
end

% Parse the different styles of input
if nargin<4 || isempty(K)
    % fitparamodel(V,model,t)
    Kpassed = false;
elseif nargin>3 && ischar(K) && ~iscell(K)
    % fitparamodel(V,model,t,'Property',Value)
    if nargin>4
        varargin = [{K} {StartParameters} varargin];
    end
    StartParameters = [];
    Kpassed = false;
elseif nargin>3 &&  ~all(size(K)>1) && ~iscell(K)
    % fitparamodel(V,model,t,StartParameters,'Property',Value)
    if nargin>4
        varargin = [{StartParameters} varargin];
    end
    StartParameters = K;
    Kpassed = false;
else
    % fitparamodel(V,model,r,K,StartParameters,'Property',Value)
    Kpassed = true;
end
isDistanceDomain = Kpassed;

% Validate input signals
if ~iscell(V)
    V = {V(:)};
end
nSignals = numel(V);
for i = 1:nSignals
    if ~iscolumn(V{i})
        V{i} = V{i}(:);
    end
    if ~isreal(V{i})
        V{i} = real(V{i});
    end
    if ~isreal(V{i})
        error('Input signal cannot be complex.')
    end
    validateattributes(V{i},{'numeric'},{'nonempty'},mfilename,'V')
end

% Validate input axis
if ~iscell(ax)
    ax = {ax(:)};
end
nAxes = numel(ax);
for i = 1:nAxes
    if ~iscolumn(ax{i})
        ax{i} = ax{i}(:);
    end
    if numel(unique(round(diff(ax{i}),6)))~=1
        error('Distance/time axis must be a monotonically increasing vector.')
    end
    % If using distance domain, only one axis must be passed
    if ~isDistanceDomain
        % Is using time-domain,control that same amount of axes as signals are passed
        if numel(V{i})~=numel(ax{i})
            error('V and t arguments must fulfill numel(t)==numel(S).')
        end
    end
end

if isDistanceDomain
    % Input #5 fitparamodel(V,model,r,K,'Property',Value)
    if nargin>4 && ischar(StartParameters)
        varargin = [{StartParameters},varargin];
        StartParameters = [];
    end
end

% Check that parametric model is passed as function handle
if ~isa(model,'function_handle')
    error('Model must be a function handle.')
end

% Get information about the parametric model
try
    % Check whether model is a DeerLab model function
    paraminfo = model();
    paraminfo = paraminfo.parameters;
    if nargin(model)==2
        model = @(ax,param,idx) model(ax,param);
    else
        model = @(ax,param,idx) model(ax,param,idx);
    end
    
catch
    % If not, require user to pass the inital values
    if ~exist('StartParameters','var') || isempty(StartParameters) || ischar(StartParameters)
        error('For this model, please provide the required inital guess parameters.')
    end
    % Wrap the function handle into a DeerLab model function
    model = paramodel(model,StartParameters,[],[],isDistanceDomain);
    paraminfo = model();
    paraminfo = paraminfo.parameters;
end

% Validate kernel
if isDistanceDomain
    if ~iscell(K)
        K = {K};
    end
    for i = 1:nSignals
        if size(K{i},1)~=numel(V{i})
            error('The number of rows in K must match the number of elements in V.')
        end
    end
    if numel(K)~=nSignals
        error('The number of kernels and signals must be equal.')
    end
end

if nargin<5 || isempty(StartParameters)
    % If user does not give parameters, use the defaults of the model
    StartParameters =  [paraminfo(:).default];
elseif nargin>4 && ischar(StartParameters)
    varargin = [{StartParameters} varargin];
    StartParameters = [paraminfo(:).default];
else
    validateattributes(StartParameters,{'numeric'},{'2d','nonempty'},mfilename,'StartParameters')
end

% Parse the optional parameters in varargin
optionalProperties = {'Solver','Algorithm','MaxIter','Verbose','MaxFunEvals',...
    'TolFun','GlobalWeights','Upper','Lower','MultiStart',...
    'ConfidenceLevel','internal::returncovariancematrix'};
[Solver,Algorithm,maxIter,Verbose,maxFunEvals,TolFun,GlobalWeights,...
    upperBounds,lowerBounds,MultiStart,ConfidenceLevel,returnCovariance] = ...
    parseoptional(optionalProperties,varargin);

% Validate optional inputs
if isempty(MultiStart)
    MultiStart = 1;
else
    validateattributes(MultiStart,{'numeric'},{'scalar','nonnegative'},mfilename,'MultiStarts')
end
if isempty(returnCovariance)
    returnCovariance = false;
else
    validateattributes(returnCovariance,{'logical'},{'nonempty'},mfilename,'returnCovariance')
end
if isempty(TolFun)
    TolFun = 1e-10;
else
    validateattributes(TolFun,{'numeric'},{'scalar','nonnegative'},mfilename,'TolFun')
end
if isempty(maxFunEvals)
    maxFunEvals = 5000;
else
    validateattributes(maxFunEvals,{'numeric'},{'scalar','nonnegative'},mfilename,'MaxFunEvals')
end
if isempty(maxIter)
    maxIter = 3000;
else
    validateattributes(maxIter,{'numeric'},{'scalar','nonnegative'},mfilename,'MaxIter')
end
if isempty(Verbose)
    Verbose = 'off';
else
    validateattributes(Verbose,{'char'},{'nonempty'},mfilename,'Verbose')
end
if isempty(ConfidenceLevel)
    ConfidenceLevel = 0.95;
else
    validateattributes(ConfidenceLevel,{'numeric'},{'nonnegative','nonempty'},mfilename,'ConfidenceLevel')
    if any(ConfidenceLevel>1 | ConfidenceLevel<0)
        error('The confidence level option must have values between 0 and 1.')
    end
end

% Solver
OptimizationToolboxInstalled = optimtoolbox_installed();
if OptimizationToolboxInstalled
    DefaultSolver = 'lsqnonlin';
    SolverList = {'lsqnonlin','nlsqbnd','lmlsqnonlin'};
else
    DefaultSolver = 'lmlsqnonlin';
    SolverList = {'nlsqbnd','lmlsqnonlin'};
end
if isempty(Solver)
    Solver = DefaultSolver;
end
validateattributes(Solver,{'char'},{'nonempty'},mfilename,'Solver')
validatestring(Solver,SolverList);
if ~ispc && strcmp(Solver,'nlsqbnd')
    error('The ''nlsqbnd'' solver is only available for Windows systems.')
end

% Algorithm
if isempty(Algorithm)
    if strcmp(Solver,'lsqnonlin')
        Algorithm = 'trust-region-reflective';
    else
        Algorithm = 'interior-point';
    end
else
    validInputs = {'levenberg-marquardt','interior-point','trust-region-reflective','active-set','sqp'};
    Algorithm = validatestring(Algorithm,validInputs);
end

% Global weights
if ~isempty(GlobalWeights)
    validateattributes(GlobalWeights,{'numeric'},{'nonnegative'})
    if numel(GlobalWeights)~=nSignals
        error('The number of global fit weights and signals must be equal.')
    end
    GlobalWeights = GlobalWeights/sum(GlobalWeights);
end


% Preparation of objective functions, parameter ranges, etc
%-------------------------------------------------------------------------------

% Parse errors in the model function, and reformat them
% model = @(ax,Parameters,idx)errorhandler(model,'modelfcn',ax,Parameters,idx);

% Get weights of different signals for global fitting
if isempty(GlobalWeights)
    Weights = globalweights(V);
else
    Weights = GlobalWeights;
end
Weights = Weights(:).';

Labels = num2cell(1:nSignals);

% Prepare upper/lower bounds on parameter search range
Ranges =  [paraminfo(:).range];
if isempty(lowerBounds)
    lowerBounds = Ranges(1:2:end-1);
end
if isempty(upperBounds)
    upperBounds = Ranges(2:2:end);
end
if any(upperBounds==realmax) || any(lowerBounds==-realmax)
    warning('Some model parameters are unbounded. Use ''Lower'' and ''Upper'' options to pass parameter boundaries.')
end
if  numel(StartParameters)~=numel(upperBounds) || ...
        numel(StartParameters)~=numel(lowerBounds)
    error('The inital guess and upper/lower boundaries must have equal length.')
end
if any(upperBounds<lowerBounds)
    error('Lower bound values cannot be larger than upper bound values.')
end

% Preprare multiple start global optimization if requested
MultiStartParameters = multistarts(MultiStart,StartParameters,lowerBounds,upperBounds);



% Configure solver
%-------------------------------------------------------------------------------
switch Solver
    
    case 'lsqnonlin'
        solverOpts = optimoptions(@lsqnonlin,'Algorithm',Algorithm,...
            'Display',Verbose,'MaxIter',maxIter,'MaxFunEvals',maxFunEvals,...
            'TolFun',TolFun,'DiffMinChange',0,'DiffMaxChange',Inf);
        solverFcn = @lsqnonlin;
        
    case 'lmlsqnonlin'
        
        solverOpts = struct('Display',Verbose,'MaxIter',maxIter,'MaxFunEvals',...
            maxFunEvals,'TolFun',TolFun);
        solverFcn = @lmlsqnonlin;
        
    case 'nlsqbnd'
        
        solverOpts = optimset('Algorithm',Algorithm,'Display',Verbose,...
            'MaxIter',maxIter,'MaxFunEvals',maxFunEvals,...
            'TolFun',TolFun,'TolCon',1e-20,...
            'DiffMinChange',1e-8,'DiffMaxChange',0.1);
        solverFcn = @nlsqbnd;
end

% Run least-squares fitting
%-------------------------------------------------------------------------------
% Disable ill-conditioned matrix warnings
warning('off','MATLAB:nearlySingularMatrix')

fvals = zeros(1,MultiStart);
parfits = cell(1,MultiStart);
for runIdx = 1:MultiStart
    
    StartParameters = MultiStartParameters(runIdx,:);
    
    [parfit,fval,~,exitflag]  = solverFcn(@ResidualsFcn,StartParameters,lowerBounds,upperBounds,solverOpts);
    if exitflag==0
        % if maxIter exceeded, doube iterations and continue
        solverOpts.MaxIter = 2*maxIter;
        solverOpts.MaxFunEvals = 2*maxFunEvals;
        [parfit,fval]  = solverFcn(@ResidualsFcn,parfit,lowerBounds,upperBounds,solverOpts);
    end
    
    fvals(runIdx) = fval;
    parfits{runIdx} = parfit;
    
end

% Find global minimum from multiple runs
[~,globmin] = min(fvals);
parfit = parfits{globmin};


% Calculate parameter confidence intervals
%-------------------------------------------------------------------------------
calcParamUncertainty = nargout>2;
if calcParamUncertainty
    
    % Compute residual vector and estimate variance from that
    residuals = ResidualsFcn(parfit);
    sigma2 = var(residuals);
    
    % Calculate numerical estimates of the Jacobian and Hessian
    % of the negative log-likelihood
    jacobian = jacobianest(@ResidualsFcn,parfit);
    hessian = jacobian.'*jacobian;
    
    % Estimate the covariance matrix by means of the inverse of Fisher information matrix
    lastwarn('');
    covmatrix = sigma2.*inv(hessian);
    % Detect if there was a 'nearly singular' warning
    [~, warnId] = lastwarn;
    if strcmp(warnId,'MATLAB:nearlySingularMatrix') || strcmp(warnId,'MATLAB:singularMatrix')
        covmatrix = sigma2.*sparse(pinv(full(hessian)));
        lastwarn('');
    end
    
    % Set significance level for confidence intervals
    alpha = 1 - ConfidenceLevel;
    p = 1 - alpha/2; % percentile
    N = numel(residuals) - numel(parfit); % degrees of freedom
    
    parci = cell(numel(p),1);
    z = zeros(numel(p),1);
    %Get the CI at requested confidence levels
    for j=1:numel(p)
        
        % Get Student's t critical value
        z(j) = t_inv(p(j),N);
        % Compute bounds of confidence intervals
        parci{j} = parfit.' + z(j)*sqrt(diag(covmatrix)).*[-1 +1];
        parci{j} = max(parci{j},lowerBounds.');
        parci{j} = min(parci{j},upperBounds.');
    end
    
    %Do not return a cell if only one confidence level is requested
    if numel(p)==1
        parci = parci{1};
    end
    
    % If wrapper functions internally request the covariance matrix, pack it up
    if returnCovariance
        parci = {parci, covmatrix, z};
    end
    
end

% Evaluate fitted model and model confidence intervals
%-------------------------------------------------------------------------------
computeFittedModel = nargout>1;
if computeFittedModel
    modelfit = cell(nAxes,1);
    for i = 1:nAxes
        modelfit{i} = model(ax{i},parfit,Labels{i});
    end
end

computeModelCI = nargout>3;
if computeModelCI
    modelci = cell(nAxes,1);
    for i = 1:nAxes
        % Compute Jacobian for time/distance-model
        jacobian = jacobianest(@(par)model(ax{i},par,Labels{i}),parfit);
        modelvariance = arrayfun(@(idx)full(jacobian(idx,:))*covmatrix*full(jacobian(idx,:)).',1:numel(ax{i})).';
        upper = modelfit{i} + z*sqrt(modelvariance);
        lower = modelfit{i} - z*sqrt(modelvariance);
        modelci{i} = [lower(:) upper(:)];
    end    
end

% Calculate goodness of fit
%-------------------------------------------------------------------------------
computeStats = nargout>4;
if computeStats
    stats = cell(nSignals,1);
    for i = 1:nSignals
        if isDistanceDomain
            Vfit = K{i}*modelfit{1};
        else
            Vfit = modelfit{i};
        end
        Ndof = numel(V{i}) - numel(StartParameters);
        stats{i} = gof(V{i},Vfit,Ndof);
    end
    
    if nSignals==1
        stats = stats{1};
    end
end

if nAxes==1
    if computeFittedModel
        modelfit = modelfit{1};
    end
    if computeModelCI
        modelci = modelci{1};
    end
end

% Set the warnings back on
warning('on','MATLAB:nearlySingularMatrix')
    
    % Function that provides vector of residuals, which is the objective
    % function for the least-squares solvers
    function r = ResidualsFcn(p)
        r_ = cell(nSignals,1);
        t = ax{1};
        for iSignal = 1:nSignals
            if nAxes>1
                t = ax{iSignal};
            end
            if isDistanceDomain
                r_{iSignal} = Weights(iSignal)*(V{iSignal}-K{iSignal}*model(t,p,iSignal));
            else
                r_{iSignal} = Weights(iSignal)*(V{iSignal}-model(t,p,iSignal));
            end
        end
        r = vertcat(r_{:});
    end
    
end
