function [Background,FitParam] = fitbackground(FitData,TimeAxis,FitTimeAxis,BckgModel,ModelParam)

if nargin<3
    error('Not enough input arguments.')
end

if nargin<5
    ModelParam = [];
end

if nargin<4 || isempty(BckgModel)
    BckgModel = 'exponential';
elseif isa(BckgModel,'function_handle')
    CustomFitModel = BckgModel;
    BckgModel = 'custom';
elseif ~isa(BckgModel,'function_handle') && ~isa(BckgModel,'char')
    error('BckgModel must be a valid ''function_handle'' or ''char'' class variable.')
else
    validateattributes(BckgModel,{'char'},{'nonempty'},mfilename,'BckgModel')
    allowedInput = {'fractal','exponential','polyexp','polynomial'};
    validatestring(BckgModel,allowedInput);
end
if iscolumn(TimeAxis)
    TimeAxis = TimeAxis';
end
if iscolumn(FitData)
    FitData = FitData';
    DataIsColumn = true;
else
    DataIsColumn = false;
end
if iscolumn(FitTimeAxis)
    FitTimeAxis = FitTimeAxis';
end
validateattributes(FitData,{'numeric'},{'2d','nonempty'},mfilename,'FitData')
validateattributes(TimeAxis,{'numeric'},{'2d','nonempty','increasing'},mfilename,'TimeAxis')
TimeAxis = abs(TimeAxis);

%--------------------------------------------------------------------------
%Memoization
%--------------------------------------------------------------------------

persistent cachedData
if isempty(cachedData)
    cachedData =  java.util.LinkedHashMap;
end
hashKey = datahash({FitData,TimeAxis,FitTimeAxis,BckgModel,ModelParam});
if cachedData.containsKey(hashKey)
    Output = cachedData.get(hashKey);
    [Background,FitParam] = java2mat(Output);
    Background = Background';
    return
end

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

%Set options for fminsearch solver
solveropts = optimset('DIsplay','off','MaxIter',5000,'MaxFunEval',10000,'TolFun',1e-15,'TolX',1e-15);

%Start with a linear fit of log-data
CostFcn = @(param)(1/2*norm(polynomial(FitTimeAxis,param) - log(FitData))^2);
PolynomialFit = fminsearch(CostFcn,[0 0],solveropts);
LinearLogFit=[-PolynomialFit(2) 1];

switch BckgModel
    
    case 'fractal'
        %Use linear log-fit as initial gues
        LinearLogFit = [-PolynomialFit(2) 1 3];
        %Fit amplitude, decay rate and fractal dimension of stretched exp.
        fminResults = fminsearch(@minimizeStretchExp,LinearLogFit,[],FitTimeAxis,FitData);
        %Limit resolution to fourth digit to avoid precission errors
        fminResults = round(fminResults,4);
        %Compute stretched exponential background
        Background = fminResults(2)*exp(-abs(fminResults(1)*TimeAxis).^(fminResults(3)/3));
        FitParam = fminResults([1,3]);
        
    case 'exponential'
        %Compute exponential background from linear log-fit
        Background = exp(polynomial(abs(TimeAxis),PolynomialFit));
        FitParam = LinearLogFit(1);
        
    case 'polynomial'
        %Fit polynomial function
        CostFcn = @(param)(1/2*norm(polynomial(FitTimeAxis,param) - FitData)^2);
        PolynomialFit = fminsearch(CostFcn,zeros(ModelParam+1,1),solveropts);
        %Compute polynomial background
        Background = polynomial(abs(TimeAxis),PolynomialFit);
        FitParam = PolynomialFit;
        
    case 'polyexp'
        %Fit polynomial function to log-data
        CostFcn = @(param)(1/2*norm(polynomial(FitTimeAxis,param) - log(FitData))^2);
        PolynomialFit = fminsearch(CostFcn,[0 0],solveropts);
        %Compute polynomial background
        Background = exp(polynomial(abs(TimeAxis),PolynomialFit));
        FitParam = PolynomialFit;
    case 'custom'
        CostFcn = @(param)(1/2*norm(CustomFitModel(FitTimeAxis,param) - FitData)^2);
        %Fit using the user-provided model
        FitParam = fminsearch(CostFcn,zeros(ModelParam,1),solveropts);
        %Compute fitted background
        Background = CustomFitModel(abs(TimeAxis),FitParam);
end

%Ensure data is real
Background=real(Background);
if DataIsColumn
    Background = Background';
end

%Store output result in the cache
Output = {Background,FitParam};
cachedData = addcache(cachedData,hashKey,Output);

%Cost function for fractal stretched exponential fit
function RMSD = minimizeStretchExp(StretchedExpParam,TimeAxis,Obervation)
if StretchedExpParam(1)<0
    RMSD=1.0e10;
    return;
end
Prediction = StretchedExpParam(2)*exp(-abs(StretchedExpParam(1)*TimeAxis).^(StretchedExpParam(3)/3));
RMSD = 1/2*norm(Prediction - Obervation)^2;
end

%Polynomial function (substitute for polyval)
function y = polynomial(x,param)
y = 0;
for i = length(param):-1:1
    y = y + param(i)*x.^(i-1);
end
end


end
