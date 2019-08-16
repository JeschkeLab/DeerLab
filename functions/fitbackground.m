%
% FITBACKGROUND Fit the background function in a signal
%
%   [B,PARAM] = FITBACKGROUND(DATA,T,TFIT,'MODEL')
%   Fits the the paramters PARAM of the N-point background function B. This
%   is done by fitting the M-point DATA on a M-point axis TFIT using a
%   model given by the string MODEL. The background is then extrapolated to
%   the N-point axis T.
%
%   [B,PARAM] = FITBACKGROUND(DATA,T,TFIT,'polynomial',ORDER)
%   For polynomial function fitting, the order of the polynomial can be
%   passed as an additional input argument.
%
%   [B,PARAM] = FITBACKGROUND(DATA,T,TFIT,@MODEL)
%   User-defined models can be fittid by passing a function handle instead
%   of a model name. 
%
% The pre-defined models in FITBACKGROUND defined by the MODEL string
% argument are the following:
%
%       'exponential' - exponential function where the decay rate if fitted
%
%       'polyexp' -  exponetial function by fitting the a linear function
%                    to the log of the signal
%
%       'fractal' - stretched exponential by fitting the decay rate and
%                   fractal dimension on the log of the signal
%
%       'polynomial' - polynomial function of order as given as an input
%
% To pass user-defined models, the MODEL argument must be a function handle
% to a function accepting two input arguments as follows
%       function myModel(AXIS,PARAM), ..., end
% where PARAM is an array containing the parameter of the model. Some of
% these models include the strexp(), sumstrexp(), prodstrexp() models
% distibuted in DeerAnalysis2.
%
%
% Copyright(C) 2019  Luis Fabregas, DeerAnalysis2
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License 3.0 as published by
% the Free Software Foundation.

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
    BckgModel = validatestring(BckgModel,allowedInput);
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
        solveropts=optimoptions(@lsqnonlin,'Algorithm','trust-region-reflective','Display','off',...
            'MaxIter',8000,'MaxFunEvals',8000,...
            'TolFun',1e-10,'DiffMinChange',1e-8,'DiffMaxChange',0.1);
        CostFcn = @(param) (sqrt(0.5)*(CustomFitModel(FitTimeAxis,param) - FitData'));
        %Get information about the parametric model
        Info = CustomFitModel();
        StartParameters =  [Info.parameters(:).default];
        Ranges =  [Info.parameters(:).range];
        LowerBounds = Ranges(1:2:end-1);
        UpperBounds = Ranges(2:2:end);
        FitParam = lsqnonlin(CostFcn,StartParameters,LowerBounds,UpperBounds,solveropts);
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
