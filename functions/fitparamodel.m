% 
% FITPARAMODEL Fits a distance distribution to one (or several) signals
%              by fitting of a parametric model.
%
%   [P,param] = FITPARAMODEL(S,K,r,@model)
%   Fitting of the N-point signal (S) to a M-point distance distribution 
%   (P) given a M-point distance axis R and NxM point kernel (K). The fitted
%   distribution corresponds to a parametric model calculated by the passed
%   function handle (@model). The fitted parameters (param) are returned as a
%   second output argument.
%
%   P = FITPARAMODEL({S1,S2,...},{K1,K2,...},r,'MODEL')
%   Passing multiple signals/kernels enables global fitting of the
%   to a single parametric model distance distribution. The global fit weights
%   are automatically computed according to their contribution to ill-posedness.
%
%   P = FITPARAMODEL(...,'Property',Values)
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



function [Distribution,FitParameters] = fitparamodel(S,K,r,Model,StartParameters,varargin)

%--------------------------------------------------------------------------
% Input Parsening & Validation
%--------------------------------------------------------------------------

%Get information about the parametric model
Info = Model();
if nargin<5 || isempty(StartParameters)
    %IF user does not give parameters, use the defaults of the model
    StartParameters =  [Info.parameters(:).default];
else
    %If user passes them, check that the number of parameters matches the model
    if length(StartParameters)~=Info.nParam
        error('The number of input parameters does not match the number of model parameters.')
    end
    validateattributes(StartParameters,{'numeric'},{'2d','nonempty'},mfilename,'StartParameters')
end

%Get optional parameters
[Solver,Algorithm,MaxIter,MaxFunEvals,TolFun,CostModel,GlobalWeights] = parseoptional({'Solver','Algorithm','MaxIter','MaxFunEvals','TolFun','CostModel','GlobalWeights'},varargin);

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
if ~iscolumn(r)
    r = r.';
end
if numel(unique(round(diff(r),12)))~=1
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
    validateattributes(S{i},{'numeric'},{'nonempty'},mfilename,'S')
end


%--------------------------------------------------------------------------
%Memoization
%--------------------------------------------------------------------------

persistent cachedData
if isempty(cachedData)
    cachedData =  java.util.LinkedHashMap;
end
hashKey = datahash({S,K,r,func2str(Model),StartParameters,varargin});
if cachedData.containsKey(hashKey)
    Output = cachedData.get(hashKey);
    [Distribution,FitParameters] = java2mat(Output);
    FitParameters = FitParameters';
    return
end

%--------------------------------------------------------------------------
% Execution
%--------------------------------------------------------------------------

%Define the cost functional of a single signal
switch CostModel
    case 'lsq'
        ModelCost = @(Parameters,K,S) (norm(K*Model(r,Parameters) - S)^2);
    case 'chisquare'
        nParam = length(StartParameters);
        ModelCost = @(Parameters,K,S) (1/(length(S) - nParam)/(noiselevel(S)^2)*sum((K*Model(r,Parameters) - S).^2));
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
LowerBounds = Ranges(1:2:end-1);
UpperBounds = Ranges(2:2:end);

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
        ModelCost = @(Parameters) (sqrt(0.5)*(K{1}*Model(r,Parameters) - S{1}));
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

%Compute fitted distance distribution
Distribution = Model(r,FitParameters);

%Store output result in the cache
Output = {Distribution,FitParameters};
cachedData = addcache(cachedData,hashKey,Output);

return

