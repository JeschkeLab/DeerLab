%
% FITPARAMODEL Fits a time- or distance-domain parametric model to one (or several) signals
%
%   [param,Vfit,paramci,fitci] = FITPARAMODEL(V,@model,t)
%   [param,Vfit,paramci,fitci] = FITPARAMODEL(V,@model,r,K)
%   Fits the N-point signal (V) to a M-point parametric model (@model) given an
%   M-point distance/time axis (r/t). For distance-domain fitting, provide
%   the NxM point kernel matrix (K). The fitted model corresponds to a parametric model
%   calculated by the passed function handle (@model). The fitted parameters (param)
%   are returned as the first output argument, their 99% confidence intervals (paramci) are
%   returned as the third output, the fitted model as the second output, and the
%   corresponding 99% confidence bands (fitci) as the fourth output.
%
%   [param,Vfit,paramci,fitci] = FITPARAMODEL(V,@model,t,param0)
%   [param,Vfit,paramci,fitci] = FITPARAMODEL(V,@model,r,K,param0)
%   The initial guess of the model parameters can be passed as a last
%   argument (param0). This is optional for DeerLab model functions. If (@model)
%   is a user-defined function handle, (param0) is required.
%
%   [param,Vfit,paramci,fitci] = FITPARAMODEL({V1,V2,___},@model,{t1,t2,___},param0)
%   [param,Vfit,paramci,fitci] = FITPARAMODEL({V1,V2,___},@model,r,{K1,K2,___},param0)
%   Pass multiple signals/kernels to enable global fitting of a single parametric
%   model to all data. The global fit weights are automatically computed according
%   to their contribution to ill-posedness.
%
%   [param,Vfit,paramci,fitci] = FITPARAMODEL(___,'Property',Values)
%   Additional (optional) arguments can be passed as name-value pairs.
%
% The properties to be passed as options can be set in any order.
%
%   'Solver' - Solver to be used to solve the minimization problems
%           'lsqnonlin' - Non-linear constrained least-squares (toolbox)
%             'fmincon' - Non-linear constrained minimization (toolbox)
%             'nlsqbnd' - Non-linear constrained least-squares (free)
%       'fminsearchbnd' - Non-linear constrained minimization (free)
%          'fminsearch' - Unconstrained minimization
%   'CostModel' - Type of fitting cost functional to use.
%                      'lsq' - Least-squares fitting
%                      'chisquared' - Chi-squared fitting (as in GLADD or DD)
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
%   'ConfidenceLevel' - Level for confidence intervals
%   'Verbose' - Display options for the solvers:
%                 'off' - no information displayed
%                 'final' - display solver exit message
%                 'iter-detailed' - display state of solver at each iteration
%               See MATLAB doc optimoptions for detailed explanation
%

% This file is a part of DeerLab. License is MIT (see LICENSE.md).
% Copyright(c) 2019-2020: Luis Fabregas, Stefan Stoll and other contributors.

function [parfit,modelfit,parci,modelci] = fitparamodel(V,model,ax,K,StartParameters,varargin)

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
if ~Kpassed
    % Check if global fitting is in use
    if iscell(V)
        K = cell(size(V));
        for i=1:length(V)
            if iscell(ax)
                K{i} = eye(length(V{i}),length(ax{i}));
            else
                K{i} = eye(length(V{i}),length(ax));
            end
        end
    else
        K = eye(length(V),length(ax));
    end
    isDistanceDomain = false;
else
    % Input #5 fitparamodel(V,model,r,K,'Property',Value)
    if nargin>4 && ischar(StartParameters)
        varargin = [{StartParameters},varargin];
        StartParameters = [];
    end
    isDistanceDomain = true;
end

% Check that parametric model is passed as function handle
if ~isa(model,'function_handle')
    error('Model must be a function handle.')
end

% Get information about the parametric model
try
    % Check whether model is a DeerLab model function
    Info = model();
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
    Info = model();
end

if nargin<5 || isempty(StartParameters)
    % If user does not give parameters, use the defaults of the model
    StartParameters =  [Info.parameters(:).default];
elseif nargin>4 && ischar(StartParameters)
    varargin = [{StartParameters} varargin];
    StartParameters = [Info.parameters(:).default];
else
    validateattributes(StartParameters,{'numeric'},{'2d','nonempty'},mfilename,'StartParameters')
end

% Parse the optional parameters in varargin
optionalProperties = {'Solver','Algorithm','MaxIter','Verbose','MaxFunEvals',...
  'TolFun','CostModel','GlobalWeights','Upper','Lower','MultiStart',...
  'ConfidenceLevel','internal::returncovariancematrix'};
[Solver,Algorithm,maxIter,Verbose,maxFunEvals,TolFun,CostModel,GlobalWeights,...
  upperBounds,lowerBounds,MultiStart,ConfidenceLevel,returnCovariance] = ...
  parseoptional(optionalProperties,varargin);

% Validate optional inputs
if isempty(CostModel)
    CostModel = 'lsq';
else
    validInputs = {'lsq','chisquare'};
    CostModel = validatestring(CostModel,validInputs);
end
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
    ConfidenceLevel = 0.99;
else
    validateattributes(ConfidenceLevel,{'numeric'},{'scalar','nonnegative','nonempty'},mfilename,'ConfidenceLevel')
    if ConfidenceLevel>1 || ConfidenceLevel<0
        error('The confidence level option must be a value between 0 and 1.')
    end
end

% Solver
OptimizationToolboxInstalled = license('test','optimization_toolbox');
if isempty(Solver) && ~OptimizationToolboxInstalled
    Solver = 'fminsearchcon';
elseif isempty(Solver) && OptimizationToolboxInstalled
    Solver = 'lsqnonlin';
else
    validateattributes(Solver,{'char'},{'nonempty'},mfilename,'Solver')
    if OptimizationToolboxInstalled
        SolverList = {'fminsearchcon','lsqnonlin','fmincon','fminsearch','nlsqbnd'};
    else
        SolverList = {'fminsearchcon','fminsearch','nlsqbnd'};
    end
    validatestring(Solver,SolverList);
end

if ~ispc && strcmp(Solver,'nlsqbnd')
    error('The ''nlsqbnd'' solver is only available for Windows systems.')
end

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

% Validate input signal and kernel
if ~iscell(V)
    V = {V(:)};
end
if ~iscell(K)
    K = {K};
end
if numel(K)~=numel(V)
    error('The number of kernels and signals must be equal.')
end
if ~isempty(GlobalWeights)
    validateattributes(GlobalWeights,{'numeric'},{'nonnegative'})
    if length(GlobalWeights)~=length(V)
        error('The number of global fit weights and signals must be equal.')
    end
    % Normalize weights
    GlobalWeights = GlobalWeights/sum(GlobalWeights);
end

for i = 1:length(V)
    if ~iscolumn(V{i})
        V{i} = V{i}.';
    end
    if ~isreal(V{i})
        V{i} = real(V{i});
    end
    if length(V{i})~=size(K{i},1)
        error('The number of rows in K must match the number of elements in V.')
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
for i = 1:length(ax)
    if ~iscolumn(ax{i})
        ax{i} = ax{i}.';
    end
    if numel(unique(round(diff(ax{i}),6)))~=1
        error('Distance/Time axis must be a monotonically increasing vector.')
    end
    % If using distance domain, only one axis must be passed
    if ~isDistanceDomain
        % Is using time-domain,control that same amount of axes as signals are passed
        if length(V{i})~=length(ax{i})
            error('V and t arguments must fulfill length(t)==length(S).')
        end
    end
    
end

%--------------------------------------------------------------------------
% Execution
%--------------------------------------------------------------------------

% Parse errors in the model function, and reformat them
% model = @(ax,Parameters,idx)errorhandler(model,'modelfcn',ax,Parameters,idx);

% Define the cost functional of a single signal
switch CostModel
    case 'lsq'
        ModelCost = @(Parameters,K,S,ax,idx) norm(K*model(ax,Parameters,idx) - S)^2;
    case 'chisquare'
        nParam = length(StartParameters);
        ModelCost = @(Parameters,K,S,ax,idx) 1/(length(S) - nParam)/(noiselevel(S)^2)*sum((K*model(ax,Parameters,idx) - S).^2);
end

% Get weights of different signals for global fitting
if isempty(GlobalWeights)
    Weights = globalweights(V);
else
    Weights = GlobalWeights;
end
Weights = Weights(:).';

Labels = num2cell(1:numel(V));
catvec = @(x) cat(1,x{:});

% Create a new handle which evaluates the model cost function for every signal
if length(ax)>1
    CostFcn = @(Parameters) (sum(Weights.*cellfun(@(K,V,t,idx)ModelCost(Parameters,K,V,t,idx),K,V,ax,Labels)));
    VecCostFcn = @(Parameters) catvec(cellfun(@(K,V,t,idx) Weights(idx)*(sqrt(0.5)*(K*model(t,Parameters,idx) - V)),K,V,ax,Labels,'UniformOutput',false));
else
    CostFcn = @(Parameters) (sum(Weights.*cellfun(@(K,V,idx)ModelCost(Parameters,K,V,ax{1},idx),K,V,Labels)));
    VecCostFcn = @(Parameters) catvec(cellfun(@(K,V,idx) Weights(idx)*(sqrt(0.5)*(K*model(ax{1},Parameters,idx) - V)),K,V,Labels,'UniformOutput',false));
end

% Prepare upper/lower bounds on parameter search
Ranges =  [Info.parameters(:).range];
if isempty(lowerBounds)
    lowerBounds = Ranges(1:2:end-1);
end
if isempty(upperBounds)
    upperBounds = Ranges(2:2:end);
end
if any(upperBounds==realmax) || any(lowerBounds==-realmax)
    warning('Some model parameters are unbounded. Use ''Lower'' and ''Upper'' options to pass parameter boundaries.')
end
if  length(StartParameters)~=length(upperBounds) || ...
    length(StartParameters)~=length(lowerBounds)
    error('The Inital guess and upper/lower boundaries must have equal length.')
end
if any(upperBounds<lowerBounds)
    error('Lower bound values cannot be larger than upper bound values.')
end

% Disable ill-conditioned matrix warnings
warning('off','MATLAB:nearlySingularMatrix')

% Preprare multiple start global optimization if requested
MultiStartParameters = multistarts(MultiStart,StartParameters,lowerBounds,upperBounds);
fvals = zeros(1,MultiStart);
fits = cell(1,MultiStart);
jacobian = [];

for runIdx = 1:MultiStart
    
    StartParameters = MultiStartParameters(runIdx,:);
    
    % Fit the parametric model...
    switch Solver
        case 'fminsearchcon'
            % ...under constraints for the parameter values range
            solverOpts=optimset('Algorithm',Algorithm,'Display',Verbose,...
                'MaxIter',maxIter,'MaxFunEvals',maxFunEvals,...
                'TolFun',TolFun,'TolCon',1e-20,...
                'DiffMinChange',1e-8,'DiffMaxChange',0.1);
            [parfit,fval,exitflag] = fminsearchcon(CostFcn,StartParameters,lowerBounds,upperBounds,[],[],[],solverOpts);
            % Check how optimization exited...
            if exitflag == 0
                % ... if maxIter exceeded (flag =0) then doube iterations and continue from where it stopped
                solverOpts=optimset('Algorithm',Algorithm,'Display',Verbose,...
                    'MaxIter',2*maxIter,'MaxFunEvals',2*maxFunEvals,...
                    'TolFun',TolFun,'TolCon',1e-10,...
                    'DiffMinChange',1e-8,'DiffMaxChange',0.1);
                [parfit,fval] = fminsearchcon(CostFcn,parfit,lowerBounds,upperBounds,[],[],[],solverOpts);
            end
            
        case 'fmincon'
            % ...under constraints for the parameter values range
            solverOpts = optimoptions(@fmincon,'Algorithm',Algorithm,'Display',Verbose,...
                'MaxIter',maxIter,'MaxFunEvals',maxFunEvals,...
                'TolFun',TolFun,'TolCon',1e-20,'StepTolerance',1e-20,...
                'DiffMinChange',1e-8,'DiffMaxChange',0.1);
            [parfit,fval,exitflag]  = fmincon(CostFcn,StartParameters,[],[],[],[],lowerBounds,upperBounds,[],solverOpts);
            % Check how optimization exited...
            if exitflag == 0
                % ... if maxIter exceeded (flag =0) then doube iterations and continue from where it stopped
                solverOpts = optimoptions(solverOpts,'MaxIter',2*maxIter,'MaxFunEvals',2*maxFunEvals,'Display',Verbose);
                [parfit,fval]  = fmincon(CostFcn,parfit,[],[],[],[],lowerBounds,upperBounds,[],solverOpts);
            end
            
        case 'lsqnonlin'
            solverOpts = optimoptions(@lsqnonlin,'Algorithm',Algorithm,'Display',Verbose,...
                'MaxIter',maxIter,'MaxFunEvals',maxFunEvals,...
                'TolFun',TolFun,'DiffMinChange',0,'DiffMaxChange',Inf);
            [parfit,fval,~,exitflag,~,~,jacobian]  = lsqnonlin(VecCostFcn,StartParameters,lowerBounds,upperBounds,solverOpts);
            
            if exitflag == 0
                % ... if maxIter exceeded (flag =0) then doube iterations and continue from where it stopped
                solverOpts = optimoptions(solverOpts,'MaxIter',2*maxIter,'MaxFunEvals',2*maxFunEvals,'Display',Verbose);
                [parfit,fval,~,~,~,~,jacobian]  = lsqnonlin(VecCostFcn,parfit,lowerBounds,upperBounds,solverOpts);
            end
            
        case 'nlsqbnd'
            solverOpts = optimset('Algorithm',Algorithm,'Display',Verbose,...
                'MaxIter',maxIter,'MaxFunEvals',maxFunEvals,...
                'TolFun',TolFun,'TolCon',1e-20,...
                'DiffMinChange',1e-8,'DiffMaxChange',0.1);
            [parfit,fval,~,exitflag] = nlsqbnd(VecCostFcn,StartParameters,lowerBounds,upperBounds,solverOpts);
            % nlsqbnd returns a column, transpose to adapt to row-style of MATLAB solvers
            parfit = parfit.';
            if exitflag == 0
                % ... if maxIter exceeded (flag =0) then doube iterations and continue from where it stopped
                solverOpts = optimset('Algorithm',Algorithm,'Display',Verbose,...
                    'MaxIter',2*maxIter,'MaxFunEvals',2*maxFunEvals,...
                    'TolFun',TolFun,'TolCon',1e-20,...
                    'DiffMinChange',1e-8,'DiffMaxChange',0.1);
                [parfit,fval] = nlsqbnd(VecCostFcn,parfit,lowerBounds,upperBounds,solverOpts);
            end
            
        case 'fminsearch'
            % ...unconstrained with all possible values
            if strcmp(Verbose,'iter-detailed')
                Verbose = 'iter';
            end
            solverOpts=optimset('Algorithm',Algorithm,'Display',Verbose,...
                'MaxIter',maxIter,'MaxFunEvals',maxFunEvals,...
                'TolFun',TolFun,'TolCon',1e-10,...
                'DiffMinChange',1e-8,'DiffMaxChange',0.1);
            [parfit,fval]  = fminsearch(CostFcn,StartParameters,solverOpts);
    end
    
    fvals(runIdx) = fval;
    fits{runIdx} = parfit;
    
end

% Find global minimum from multiple runs
[~,globmin] = min(fvals);
parfit = fits{globmin};

if nargout>2
    
    % Numerically estimate the Jacobian if not done by MATLAB's lsqnonlin
    jacobian = jacobianest(VecCostFcn,parfit);
    hessian = jacobian'*jacobian;
    % Compute residual vector
    residual = VecCostFcn(parfit);
    lastwarn('');
    % Set significance level for confidence intervals
    alpha = 1 - ConfidenceLevel;
    %Get Student's T critical value
    critical = tinv(1 - alpha/2,length(residual) - numel(parfit));
    % Estimate the covariance matrix by means of the inverse of Fisher information matrix
    covmatrix = var(residual).*inv(hessian);
    % Detect if there was a 'nearly singular' warning
    [~, warnId] = lastwarn;
    if strcmp(warnId,'MATLAB:nearlySingularMatrix') || strcmp(warnId,'MATLAB:singularMatrix')
        covmatrix = var(residual).*sparse(pinv(full(hessian)));
        lastwarn('');
    end
    % Compute upper/lower confidence intervals
    parci = nan(numel(parfit),2);
    parci(:,1) = parfit - critical*sqrt(diag(covmatrix).');
    parci(:,2) = parfit + critical*sqrt(diag(covmatrix).');
    
    %If wrapper functions internally request the covariance matrix, pack it up
    if returnCovariance
        tmp{1} = parci;
        tmp{2} = covmatrix;
        tmp{3} = critical;
        parci = tmp;
    end
    
end

% Compute fitted parametric model and confidence bands if requested
if nargout>1
    modelfit = cell(numel(ax),1);
    modelci = cell(numel(ax),1);
    for i = 1:length(ax)
        modelfit{i} = model(ax{i},parfit,Labels{i});
        if nargout>3
            %Compute Jacobian for time/distance-model
            jacobian = jacobianest(@(par)model(ax{i},par,Labels{i}),parfit);
            modelvariance = arrayfun(@(idx)full(jacobian(idx,:))*covmatrix*full(jacobian(idx,:)).',1:numel(ax{i})).';
            upper = modelfit{i} + critical*sqrt(modelvariance);
            lower = modelfit{i} - critical*sqrt(modelvariance);
            modelci{i} = [upper(:) lower(:)];
        end
    end
    if length(modelfit)==1
        modelfit = modelfit{1};
        if nargout>3
            modelci = modelci{1};
        end
    end
end

% Set the warnings back on
warning('on','MATLAB:nearlySingularMatrix')

end
