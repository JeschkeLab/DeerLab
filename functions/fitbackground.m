%
% FITBACKGROUND Fit the background function in a signal
%
%   [B,lambda,param] = FITBACKGROUND(S,t,@model)
%   Fits the background (B) and the modulation depth (lambda) to a
%   time-domain signal (S) and time-axis (t) based on a given time-domain
%   parametric model (@model). The fitted parameters of the model are
%   returned as a last output argument.
%
%   [B,lambda,param] = FITBACKGROUND(S,t,@model,tstart)
%   The time at which the background starts to be fitted can be passed as a
%   an additional argument. 
%   
%   [B,lambda,param] = FITBACKGROUND(S,t,@model,[tstart tend])
%   The start and end times of the fitting can be specified by passing a
%   two-element array as the argument. If tend is not specified, the end of
%   the signal is selected as the default.
%
%   [B,lambda,param] = FITBACKGROUND(...,'Property',Value)
%   Additional (optional) arguments can be passed as property-value pairs.
%
% The properties to be passed as options can be set in any order.
%
%   'LogFit' - Specifies whether to fit the log of the signal (default: false)
%
% Copyright(C) 2019  Luis Fabregas, DeerAnalysis2
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License 3.0 as published by
% the Free Software Foundation.

function [Background,ModDepth,FitParam] = fitbackground(Data,TimeAxis,BckgModel,FitDelimiter,varargin)


[LogFit] = parseoptional({'LogFit'},varargin);


if nargin<3
    error('Not enough input arguments.')
end

if nargin<4
    FitDelimiter = minmax(TimeAxis);
elseif length(FitDelimiter) == 1
    FitDelimiter(2) = max(TimeAxis);
elseif length(FitDelimiter) > 2
    error('The 4th argument cannot exceed two elements.')
end

if FitDelimiter(2)<FitDelimiter(1)
    error('The fit start time cannot exceed the fit end time.')
end

if ~isa(BckgModel,'function_handle')
   error('The background model must be a valid function handle.') 
end

if ~iscolumn(TimeAxis)
    TimeAxis = TimeAxis.';
end

if isempty(LogFit)
   LogFit = false; 
end

if iscolumn(Data)
    DataIsColumn = true;
else
    Data = Data.';
    DataIsColumn = false;
end

validateattributes(FitDelimiter,{'numeric'},{'2d','nonempty'},mfilename,'FitDelimiter')
validateattributes(Data,{'numeric'},{'2d','nonempty'},mfilename,'Data')
validateattributes(TimeAxis,{'numeric'},{'2d','nonempty','increasing'},mfilename,'TimeAxis')

%--------------------------------------------------------------------------
%Memoization
%--------------------------------------------------------------------------

persistent cachedData
if isempty(cachedData)
    cachedData =  java.util.LinkedHashMap;
end
hashKey = datahash({Data,TimeAxis,BckgModel,FitDelimiter});
if cachedData.containsKey(hashKey)
    Output = cachedData.get(hashKey);
    [Background,ModDepth,FitParam] = java2mat(Output);
    %Java does not recognize columns
    Background = Background';
    return
end

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

%Find the position to limit fit
FitStartTime = FitDelimiter(1);
FitEndTime = FitDelimiter(2);
[~,FitStartPos] = min(abs(TimeAxis - FitStartTime)); 
[~,FitEndPos] = min(abs(TimeAxis - FitEndTime));

%Limit the time axis and the data to fit
FitTimeAxis = TimeAxis(FitStartPos:FitEndPos);
FitData = Data(FitStartPos:FitEndPos);

%Use absolute time scale to ensure proper fitting of negative-time data
FitTimeAxis = abs(FitTimeAxis);

%Prepare minimization problem solver
solveropts = optimoptions(@lsqnonlin,'Algorithm','trust-region-reflective','Display','off',...
    'MaxIter',8000,'MaxFunEvals',8000,...
    'TolFun',1e-10,'DiffMinChange',1e-8,'DiffMaxChange',0.1);

%Construct cost functional for minimization
if LogFit
    CostFcn = @(param)(sqrt(1/2)*(log((1 - param(1) + eps)*BckgModel(FitTimeAxis,param(2:end))) - log(FitData)));
else
    CostFcn = @(param)(sqrt(1/2)*((1 - param(1))*BckgModel(FitTimeAxis,param(2:end)) - FitData));
end

%Initiallize StartParameters (1st element is modulation depth)
StartParameters(1) = 0.5;
LowerBounds(1) = 0;
UpperBounds(1) = 1;

%Get information about the time-domain parametric model
Info = BckgModel();
StartParameters(2:1 + Info.nParam) =  [Info.parameters(:).default];
Ranges =  [Info.parameters(:).range];
LowerBounds(2:1 + Info.nParam) = Ranges(1:2:end-1);
UpperBounds(2:1 + Info.nParam) = Ranges(2:2:end);
FitParam = lsqnonlin(CostFcn,StartParameters,LowerBounds,UpperBounds,solveropts);

%Get the fitted modulation depth
ModDepth = FitParam(1);

%Extrapolate fitted background to whole time axis
Background = BckgModel(abs(TimeAxis),FitParam(2:end));

%Ensure data is real
Background = real(Background);
if ~DataIsColumn
    Background = Background';
end

%Store output result in the cache
Output = {Background,ModDepth,FitParam};
cachedData = addcache(cachedData,hashKey,Output);

end
