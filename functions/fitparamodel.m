function [Distribution,FitParameters] = fitparamodel(Signal,Kernel,DistanceAxis,Model,StartParameters,varargin)

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
[Solver,Algorithm,MaxIter,MaxFunEvals,TolFun,CostModel] = parseoptional({'Solver','Algorithm','MaxIter','MaxFunEvals','TolFun','CostModel'},varargin);

if isempty(CostModel)
    CostModel = 'lsq';
else
    validInputs = {'lsq','chisquare'};
    validatestring(CostModel,validInputs);
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
    validatestring(Algorithm,validInputs);
end
if ~iscolumn(DistanceAxis)
    DistanceAxis = DistanceAxis.';
end
if ~iscell(Signal)
    Signal = {Signal};
end
if ~iscell(Kernel)
    Kernel = {Kernel};
end
if length(Kernel)~=length(Signal)
    error('The number of kernels and signals must be equal.')
end
if length(Signal)>1 && strcmp(Solver,'lsqnonlin')
    Solver = 'fmincon';
    Algorithm = 'interior-point';
end
for i=1:length(Signal)
    if ~iscolumn(Signal{i})
        Signal{i} = Signal{i}.';
    end
    if ~isreal(Signal{i})
        Signal{i} = real(Signal{i});
    end
    if length(Signal{i})~=size(Kernel{i},1)
        error('Kernel and signal arguments must fulfill size(Kernel,1)==length(Signal).')
    end
    validateattributes(Signal{i},{'numeric'},{'nonempty'},mfilename,'Signal')
end


%--------------------------------------------------------------------------
%Memoization
%--------------------------------------------------------------------------

persistent cachedData
if isempty(cachedData)
    cachedData =  java.util.LinkedHashMap;
end
hashKey = datahash({Signal,Kernel,DistanceAxis,func2str(Model),StartParameters,varargin});
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
        ModelCost = @(Parameters,Kernel,Signal) (norm(Kernel*Model(DistanceAxis,Parameters) - Signal)^2);
    case 'chisquare'
        nParam = length(StartParameters);
        ModelCost = @(Parameters,Kernel,Signal) (1/(length(Signal) - nParam)/(noiselevel(Signal)^2)*sum((Kernel*Model(DistanceAxis,Parameters) - Signal).^2));
end

%Get weights of different signals for global fitting
Weights = globalweights(Signal);
%Create a new handle which evaluates the model cost function for every signal
CostFcn = @(Parameters) (sum(Weights.*cellfun(@(x,y)ModelCost(Parameters,x,y),Kernel,Signal)));

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
        ModelCost = @(Parameters) (sqrt(0.5)*(Kernel{1}*Model(DistanceAxis,Parameters) - Signal{1}));
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
Distribution = Model(DistanceAxis,FitParameters);

%Store output result in the cache
Output = {Distribution,FitParameters};
cachedData = addcache(cachedData,hashKey,Output);

return

