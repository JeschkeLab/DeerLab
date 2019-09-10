%
% FITPARAMODEL Fits a distance distribution to one (or several) signals
%              by fitting of a parametric model.
%
%   [fit,param] = FITPARAMODEL(S,@model,t)
%   [fit,param] = FITPARAMODEL(S,@model,r,K)
%   Fitting of the N-point signal (S) to a M-point parametric model
%   (fit) given a M-point distance/time axis (r/t). For distance-domain fitting
%   the NxM point kernel (K). The fitted model corresponds to a parametric model
%   calculated by the passed function handle (@model). The fitted parameters (param)
%   are returned as a second output argument.
%
%   [fit,param] = FITPARAMODEL(S,@model,t,param0)
%   [fit,param] = FITPARAMODEL(S,@model,r,K,param0)
%   The initial guess of the model parameters can be passed as a last
%   argument (param0). If (@model) is a user-defined function handle, it is
%   required to pass (param0) as an arugment.
%
%   [fit,param] = FITPARAMODEL({S1,S2,...},@model,t,param0)
%   [fit,param] = FITPARAMODEL({S1,S2,...},@model,r,{K1,K2,...},param0)
%   Passing multiple signals/kernels enables global fitting of the
%   to a single parametric model distance distribution. The global fit weights
%   are automatically computed according to their contribution to ill-posedness.
%
%   [fit,param] = FITPARAMODEL(...,'Property',Values)
%   Additional (optional) arguments can be passed as property-value pairs.
%
% The properties to be passed as options can be set in any order.
%
%   'Solver' - Solver to be used to solve the minimization problems
%                      'lsqnonlin' - Non-linear constrained least-squares
%                      'fmincon' - Non-linear constrained minimization
%                      'fminsearch' - Unconstrained minimization
%
%   'CostModel' - Type of fitting cost functional to use.
%                      'lsq' - Least-squares fitting
%                      'chisquared' - Chi-squared fitting (as in GLADD or DD)
%
%   'GlobalWeights' - Array of weighting coefficients for the individual signals in
%                     global fitting regularization.
%
%   'Algorithm' - Algorithm to be used by the solvers (see fmincon or
%                 lsqnonlin documentation)
%
%   'Lower' - Lower bound on the model parameters
%
%   'Upper' - Upper bound on the model parameters
%
%   'TolFun' - Optimizer function tolerance
%
%   'MaxIter' - Maximum number of optimizer iterations
%
%   'MaxFunEvals' - Maximum number of optimizer function evaluations
%
%
% Copyright(C) 2019  Luis Fabregas, DeerAnalysis2
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License 3.0 as published by
% the Free Software Foundation.



function [P,FitParameters] = fitparamodel(S,model,ax,K,StartParameters,varargin)

%--------------------------------------------------------------------------
% Input Parsening & Validation
%--------------------------------------------------------------------------

%Parse the different styles of input

%Input #1 fitparamodel(S,model,t)
if nargin<4 || isempty(K)
    Knotpassed = true;
%Input #2 fitparamodel(S,model,t,'Property',Value)
elseif nargin>3 && ischar(K) && ~iscell(K)
    if nargin>4
        varargin = [{K} {StartParameters} varargin];
    end
    StartParameters = [];
    Knotpassed = true;
%Input #3 fitparamodel(S,model,t,StartParameters,'Property',Value)
elseif nargin>3 &&  ~all(size(K)>1) && ~iscell(K)
    if nargin>4
        varargin = [{StartParameters} varargin];
    end
    StartParameters = K;
    Knotpassed = true;
%Input #4 fitparamodel(S,model,r,K,StartParameters,'Property',Value)
else
    Knotpassed = false;
end
if Knotpassed
    %Check if global fitting is in use
    if iscell(S)
        K = cell(size(S));
        for i=1:length(S)
         K{i} = eye(length(S{i}),length(ax));
        end
    else
        K = eye(length(S),length(ax));
    end
    isDistanceDomain = false;
else
%Input #5 fitparamodel(S,model,r,K,'Property',Value)
    if nargin>4 && ischar(StartParameters)
        varargin = [{StartParameters},varargin];
        StartParameters = [];
    end
    isDistanceDomain = true;
end

if ~isa(model,'function_handle')
   error('Model must be a valid function handle.') 
end

%Get information about the parametric model
try
    %Check whether model is a DeerAnalysis model...
    Info = model();
catch
    %... if not, then user is required to pass the inital values
    if isempty(StartParameters) || ischar(StartParameters)
        error('When using a user-defined function, the inital guess parameters are required.')
    end
    %If passed, then transform the function handle to valid parametric model
    model = paramodel(model,StartParameters,[],[],isDistanceDomain);
    Info = model();
end
if nargin<5 || isempty(StartParameters)
    %If user does not give parameters, use the defaults of the model
    StartParameters =  [Info.parameters(:).default];
elseif nargin > 4 && ischar(StartParameters)
    varargin = [{StartParameters} varargin];
    StartParameters = [Info.parameters(:).default];
else
    validateattributes(StartParameters,{'numeric'},{'2d','nonempty'},mfilename,'StartParameters')
end

%Parse the optional parameters in the varargin
[Solver,Algorithm,MaxIter,MaxFunEvals,TolFun,CostModel,GlobalWeights,UpperBounds,LowerBounds] = parseoptional(...
    {'Solver','Algorithm','MaxIter','MaxFunEvals','TolFun','CostModel','GlobalWeights','Upper','Lower'},varargin);

%Validate all inputs
if isempty(CostModel)
    CostModel = 'lsq';
else
    validInputs = {'lsq','chisquare'};
    CostModel = validatestring(CostModel,validInputs);
end

if isempty(TolFun)
    TolFun = 1e-10;
else
    validateattributes(TolFun,{'numeric'},{'scalar','nonnegative'})
end

if isempty(MaxFunEvals)
    MaxFunEvals = 5000;
else
    validateattributes(MaxFunEvals,{'numeric'},{'scalar','nonnegative'})
end

if isempty(MaxIter)
    MaxIter = 3000;
else
    validateattributes(MaxIter,{'numeric'},{'scalar','nonnegative'})
end

if isempty(Solver)
    Solver = 'lsqnonlin';
else
    validateattributes(Solver,{'char'},{'nonempty'},mfilename,'Solver')
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
if ~iscolumn(ax)
    ax = ax.';
end
if numel(unique(round(diff(ax),6)))~=1
    error('Distance axis must be a monotonically increasing vector.')
end
if ~iscell(S)
    S = {S};
end
if ~iscell(K)
    K = {K};
end
if length(K)~=length(S)
    error('The number of kernels and signals must be equal.')
end
if ~isempty(GlobalWeights)
    validateattributes(GlobalWeights,{'numeric'},{'nonnegative'})
    if length(GlobalWeights) ~= length(S)
        error('The same number of global fit weights as signals must be passed.')
    end
    if sum(GlobalWeights)~=1
        error('The sum of the global fit weights must equal 1.')
    end
end
if length(S)>1 && strcmp(Solver,'lsqnonlin')
    Solver = 'fmincon';
    Algorithm = 'interior-point';
end
for i=1:length(S)
    if ~iscolumn(S{i})
        S{i} = S{i}.';
    end
    if ~isreal(S{i})
        S{i} = real(S{i});
    end
    if length(S{i})~=size(K{i},1)
        error('K and signal arguments must fulfill size(K,1)==length(S).')
    end
    if ~isreal(S{i})
        error('Input signal cannot be complex.')
    end
    validateattributes(S{i},{'numeric'},{'nonempty'},mfilename,'S')
end

if isDistanceDomain
    %Convert distance axis to nanometers if given in Angstrom
    if ~isnanometer(ax)
        ax = ax/10;
    end 
else
    % Convert time axis to microseconds if given in nanoseconds
    usesNanoseconds = mean(diff(ax))>=0.5;
    if usesNanoseconds
        ax = round(ax)/1000; % ns->us
    end
end


%--------------------------------------------------------------------------
% Execution
%--------------------------------------------------------------------------

%Define the cost functional of a single signal
switch CostModel
    case 'lsq'
        ModelCost = @(Parameters,K,S) (norm(K*model(ax,Parameters) - S)^2);
    case 'chisquare'
        nParam = length(StartParameters);
        ModelCost = @(Parameters,K,S) (1/(length(S) - nParam)/(noiselevel(S)^2)*sum((K*model(ax,Parameters) - S).^2));
end

%Get weights of different signals for global fitting
if isempty(GlobalWeights)
    Weights = globalweights(S);
else
    Weights = GlobalWeights;
end
%Create a new handle which evaluates the model cost function for every signal
CostFcn = @(Parameters) (sum(Weights.*cellfun(@(x,y)ModelCost(Parameters,x,y),K,S)));

%Prepare upper/lower bounds on parameter search
Ranges =  [Info.parameters(:).range];
if isempty(LowerBounds)
    LowerBounds = Ranges(1:2:end-1);
end
if isempty(UpperBounds)
    UpperBounds = Ranges(2:2:end);
end
if any(UpperBounds==realmax) || any(LowerBounds == -realmax)
   warning('Some model parameters are unbounded. Use ''Lower'' and ''Upper'' options to pass parameter boundaries')
end
if any(length(StartParameters)~=length(UpperBounds) & length(StartParameters)~=length(LowerBounds))
   error('The Inital guess and upper/lower boundaries must have equal length); ') 
end
if any(UpperBounds<LowerBounds)
   error('Lower bound values cannot be larger than upper bound values.') 
end

%Disable ill-conditioned matrix warnings
warning('off','MATLAB:nearlySingularMatrix')

%Fit the parametric model...
switch Solver
    case 'fmincon'
        %...under constraints for the parameter values range
        solverOpts=optimoptions(@fmincon,'Algorithm',Algorithm,'Display','off',...
            'MaxIter',MaxIter,'MaxFunEvals',MaxFunEvals,...
            'TolFun',TolFun,'TolCon',1e-10,...
            'DiffMinChange',1e-8,'DiffMaxChange',0.1);
        [FitParameters,~,exitflag]  = fmincon(CostFcn,StartParameters,[],[],[],[],LowerBounds,UpperBounds,[],solverOpts);
        %Check how optimization exited...
        if exitflag == 0
            %... if maxIter exceeded (flag =0) then doube iterations and continue from where it stopped
            solverOpts=optimoptions(solverOpts,'MaxIter',2*MaxIter,'MaxFunEvals',2*MaxFunEvals);
            [FitParameters]  = fmincon(CostFcn,FitParameters,[],[],[],[],LowerBounds,UpperBounds,[],solverOpts);
        end
        
    case 'lsqnonlin'
        
        solverOpts=optimoptions(@lsqnonlin,'Algorithm',Algorithm,'Display','off',...
            'MaxIter',MaxIter,'MaxFunEvals',MaxFunEvals,...
            'TolFun',TolFun,'DiffMinChange',1e-8,'DiffMaxChange',0.1);
        ModelCost = @(Parameters) (sqrt(0.5)*(K{1}*model(ax,Parameters) - S{1}));
        [FitParameters,~,~,exitflag]  = lsqnonlin(ModelCost,StartParameters,LowerBounds,UpperBounds,solverOpts);
        if exitflag == 0
            %... if maxIter exceeded (flag =0) then doube iterations and continue from where it stopped
            solverOpts=optimoptions(solverOpts,'MaxIter',2*MaxIter,'MaxFunEvals',2*MaxFunEvals);
            [FitParameters]  = lsqnonlin(ModelCost,FitParameters,LowerBounds,UpperBounds,solverOpts);
        end
        
    case 'fminsearch'
        %...unconstrained with all possible values
        solverOpts=optimset('Algorithm',Algorithm,'Display','off',...
            'MaxIter',MaxIter,'MaxFunEvals',MaxFunEvals,...
            'TolFun',TolFun,'TolCon',1e-10,...
            'DiffMinChange',1e-8,'DiffMaxChange',0.1);
        FitParameters  = fminsearch(CostFcn,StartParameters,solverOpts);
end

%Set the warnings back on
warning('on','MATLAB:nearlySingularMatrix')

%Compute fitted distance distribution
P = model(ax,FitParameters);

return

