function [Distribution,FitParameters] = fitparamodel(Signal,Kernel,DistanceAxis,Model,StartParameters,varargin)

% Input Parsening & Validation
%========================================================

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
[Constrained,Algorithm,MaxIter,MaxFunEvals,TolFun] = parseoptional({'Constrained','Algorithm','MaxIter','MaxFunEvals','TolFun'},varargin);

if isempty(Algorithm)
   Algorithm = 'interior-point';
else
    validInputs = {'interior-point','trust-region-reflective','active-set','sqp'};
    validatestring(Algorithm,validInputs)
end

if isempty(TolFun)
   TolFun = 1e-10;
else
    validateattributes(TolFun,{'numeric'},{'scalar','nonnegative'})
end

if isempty(MaxFunEvals)
   MaxFunEvals = 3000;
else
    validateattributes(MaxFunEvals,{'numeric'},{'scalar','nonnegative'})
end

if isempty(MaxIter)
   MaxIter = 3000;
else
    validateattributes(MaxIter,{'numeric'},{'scalar','nonnegative'})
end

if isempty(Constrained)
    Constrained = true;
else
    validateattributes(Constrained,{'logical'},{'scalar'},mfilename,'Constrained')
end

% Execution
%========================================================

%Define the cost functional
ModelCost = @(Parameters) (1/2*norm(Kernel*Model(DistanceAxis,Parameters) - Signal)^2);

%Fit the parametric model...
if Constrained
    %...under constraints for the parameter values range
    solverOpts=optimoptions(@fmincon,'Algorithm',Algorithm,'Display','off',...
                            'MaxIter',MaxIter,'MaxFunEvals',MaxFunEvals,...
                            'TolFun',TolFun,'TolCon',1e-10,...
                            'DiffMinChange',1e-8,'DiffMaxChange',0.1);
    Ranges =  [Info.parameters(:).range];
    LowerBounds = Ranges(1:2:end-1);
    UpperBounds = Ranges(2:2:end);
    FitParameters  = fmincon(ModelCost,StartParameters,[],[],[],[],LowerBounds,UpperBounds,[],solverOpts);
else
    %...unconstrained with all possible values
    solverOpts=optimset('Algorithm',Algorithm,'Display','off',...
                            'MaxIter',MaxIter,'MaxFunEvals',MaxFunEvals,...
                            'TolFun',TolFun,'TolCon',1e-10,...
                            'DiffMinChange',1e-8,'DiffMaxChange',0.1);
    FitParameters  = fminsearch(ModelCost,StartParameters,solverOpts);
end

%Compute fitted distance distribution
Distribution = Model(DistanceAxis,FitParameters);

return

