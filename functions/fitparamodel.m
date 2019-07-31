function [Distribution,FitParameters] = fitparamodel(Signal,Kernel,DistanceAxis,Model,StartParameters,varargin)

Info = Model();
if nargin<5 || isempty(StartParameters)
    %IF user does not give parameters, use the defaults of the model
    StartParameters =  [Info.parameters(:).default];
else
    %If user passes them, check that the number of parameters matches the model
    if length(StartParameters)~=Info.nParam
        error('The number of input parameters does not match the model parameters.')
    end
    validateattributes(StartParameters,{'numeric'},{'2d','nonempty'},mfilename,'StartParameters')
end

[Constrained] = parseoptional({'Constrained'},varargin);

if isempty(Constrained)
    Constrained = true;
else
    validateattributes(Constrained,{'logical'},{'scalar'},mfilename,'Constrained')
end

%Define the cost functional
ModelCost = @(Parameters) (1/2*norm(Kernel*Model(DistanceAxis,Parameters) - Signal)^2);

%Fit the parametric model...
if Constrained
    %...under constraints for the parameter values range
    solverOpts=optimoptions(@fmincon,'Algorithm','Active-Set','Display','off',...
                            'MaxIter',300,'MaxFunEvals',3000,...
                            'TolFun',1e-10,'TolCon',1e-10,...
                            'DiffMinChange',1e-8,'DiffMaxChange',0.1);
    Ranges =  [Info.parameters(:).range];
    LowerBounds = Ranges(1:2:end-1);
    UpperBounds = Ranges(2:2:end);
    FitParameters  = fmincon(ModelCost,StartParameters,[],[],[],[],LowerBounds,UpperBounds,[],solverOpts);
else
    %...unconstrained with all possible values
    FitParameters  = fminsearch(ModelCost,StartParameters);
end

%Compute fitted distance distribution
Distribution = Model(DistanceAxis,FitParameters);

return







