%
% BACKGROUNDSTART Computes the optimal point for the background fitting start
% 
%   [t0,pos] = BACKGROUNDSTART(S,t,'model')
%   Returns the optimal start time (t0) and corresponding array index (pos)
%   at which to start fitting the background function corresponding to the
%   model given by 'model';
% 
%   [t0,pos] = BACKGROUNDSTART(...,'Property',Value)
%   Additional (optional) arguments can be passed as property-value pairs.
%
% The properties to be passed as options can be set in any order.
%
%   'RelSearchStart' - Relative position at which the background start
%                      search starts (default=0.1).
%
%   'RelSearchEnd' - Relative position at which the backgrund start search
%                    stops (default=0.6).
%
%   'EndCutOff' - Number of points to ignore from the end of the signal
%
%   'ModelParam' - Parameters for the background models
%
%   For further property-value pair options see "help fitbackground"
%
% Adapted from Gunnar Jeschke
%
% Copyright(C) 2019  Luis Fabregas, DeerAnalysis2
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License 3.0 as published by
% the Free Software Foundation.

function [FitStartTime,FitStartPos] = backgroundstart(Signal,TimeAxis,BckgModel,varargin)


%--------------------------------------------------------------------------
% Parse & Validate Required Input
%--------------------------------------------------------------------------
if nargin<2
    error('Not enough input arguments.')
end



if nargin<3 || isempty(BckgModel)
    BckgModel = 'exponential';
end

if iscolumn(TimeAxis)
    TimeAxis = TimeAxis';
end
if iscolumn(Signal)
    Signal = Signal';
end
validateattributes(Signal,{'numeric'},{'2d','nonempty'},mfilename,'FitData')
validateattributes(TimeAxis,{'numeric'},{'2d','nonempty','nonnegative','increasing'},mfilename,'TimeAxis')



%--------------------------------------------------------------------------
% Parse & Validate Optional Input
%--------------------------------------------------------------------------
%Check if user requested some options via name-value input
[RelSearchStart,RelSearchEnd,EndCutoffPos,ModelParam] = parseoptional({'RelSearchStart','RelSearchEnd','EndCutoffPos','ModelParam'},varargin);

if isempty(RelSearchStart)
    RelSearchStart = 0.1;
end

if isempty(RelSearchEnd)
    RelSearchEnd = 0.6;
end
if  isempty(EndCutoffPos)
    EndCutoffPos = length(TimeAxis);
else
    validateattributes(EndCutoffPos,{'numeric'},{'scalar','nonempty'},mfilename,'EndCutoffPos')
end

%--------------------------------------------------------------------------
%Memoization
%--------------------------------------------------------------------------

persistent cachedData
if isempty(cachedData)
    cachedData =  java.util.LinkedHashMap;
end
hashKey = datahash({Signal,TimeAxis,BckgModel,varargin});
if cachedData.containsKey(hashKey)
    Output = cachedData.get(hashKey);
    [FitStartTime,FitStartPos] = java2mat(Output);
    return
end

%--------------------------------------------------------------------------
% Adaptive background correction start search
%--------------------------------------------------------------------------
%Get zero-time position
[~,ZeroTimePosition] = min(abs(TimeAxis));
TimeAxis = TimeAxis(ZeroTimePosition:EndCutoffPos);
Signal = real(Signal((ZeroTimePosition:EndCutoffPos)));
Length = length(TimeAxis);

APTkernel = aptkernel(TimeAxis);

%Get APT kernel data
Kernel = APTkernel.Base;
NormConstant = APTkernel.NormalizationFactor;
APT_TimeAxis = APTkernel.TimeAxis;
Crosstalk = APTkernel.Crosstalk;

%Search for background fit start
StartPosMin = round(RelSearchStart*Length);
if StartPosMin<1
    StartPosMin=1;
end
StartPosMax = round(RelSearchEnd*Length);
if StartPosMax<5
    StartPosMax=5;
end

%Preallocate Merit variable
Merit = zeros(1,StartPosMax-StartPosMin);

for FitStartPos = StartPosMin:StartPosMax
    
    %Define data to be fitted according to current background start
    FitTimeAxis = TimeAxis(FitStartPos:length(TimeAxis));
    FitData = Signal(FitStartPos:length(Signal));
    
    %Fit the background with current start
    Background = fitbackground(FitData,TimeAxis,FitTimeAxis,BckgModel,ModelParam);
    
    %Correct the background from the from factor
    FormFactor = Signal - Background;
    FormFactor = FormFactor./Background;
    FormFactor = FormFactor/max(FormFactor);
    
    %Perform APT on background-corrected signal
    [FreqDimension,~] = size(Kernel);
    FreqDistribution=zeros(1,FreqDimension);
    for k=1:FreqDimension % sum in eqn [21]
        FreqDistribution(k)=FreqDistribution(k)+sum(Kernel(k,:).*FormFactor.*APT_TimeAxis)/NormConstant(k);
    end
    APTdistribution = Crosstalk\FreqDistribution';
    
    %Get merit value for this background start value
    Merit(FitStartPos - StartPosMin + 1) = sum(abs(APTdistribution(1:3)));
end

[~,OptStartPos] = min(Merit);
FitStartPos = OptStartPos + StartPosMin - 1;
FitStartTime = TimeAxis(FitStartPos);

%Store output result in the cache
Output = {FitStartTime,FitStartPos};
cachedData = addcache(cachedData,hashKey,Output);
