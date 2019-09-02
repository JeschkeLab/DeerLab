%
% BACKGROUNDSTART Computes the optimal point for the background fitting start
% 
%   [t0,pos] = BACKGROUNDSTART(V,t,@model)
%   Returns the optimal start time (t0) and corresponding array index (pos)
%   at which to start fitting a background model (@model) to the
%   data (V).
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

function [FitStartTime,FitStartPos] = backgroundstart(S,t,BckgModel,varargin)


%--------------------------------------------------------------------------
% Parse & Validate Required Input
%--------------------------------------------------------------------------
if nargin<2
    error('Not enough input arguments.')
end

if ~isa(BckgModel,'function_handle')
   error('The background model must be a valid function handle.') 
end

if iscolumn(t)
    t = t';
end
if iscolumn(S)
    S = S';
end
validateattributes(S,{'numeric'},{'2d','nonempty'},mfilename,'FitData')
validateattributes(t,{'numeric'},{'2d','nonempty','increasing'},mfilename,'t')

%Convert time step to microseconds if given in nanoseconds
usesNanoseconds = mean(diff(t))>=0.5;
if usesNanoseconds
    t = t/1000; % ns->us
end

%--------------------------------------------------------------------------
% Parse & Validate Optional Input
%--------------------------------------------------------------------------
%Check if user requested some options via name-value input
[RelSearchStart,RelSearchEnd,EndCutoffPos] = parseoptional({'RelSearchStart','RelSearchEnd','EndCutoffPos'},varargin);

if isempty(RelSearchStart)
    RelSearchStart = 0.1;
end

if isempty(RelSearchEnd)
    RelSearchEnd = 0.6;
end
if  isempty(EndCutoffPos)
    EndCutoffPos = length(t);
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
hashKey = datahash({S,t,BckgModel,varargin});
if cachedData.containsKey(hashKey)
    Output = cachedData.get(hashKey);
    [FitStartTime,FitStartPos] = java2mat(Output);
    return
end

%--------------------------------------------------------------------------
% Adaptive background correction start search
%--------------------------------------------------------------------------
%Get zero-time position
[~,ZeroTimePosition] = min(abs(t));
t = t(ZeroTimePosition:EndCutoffPos);
S = real(S((ZeroTimePosition:EndCutoffPos)));
Length = length(t);

APTkernel = aptkernel(t);

%Get APT kernel data
K = APTkernel.Base;
NormConstant = APTkernel.NormalizationFactor;
APT_t = APTkernel.t;
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
    FitStart = t(FitStartPos);
    
    %Fit the background with current start
    S = S/max(S);
    [B,ModDepth] = fitbackground(S,t,BckgModel,FitStart);

    %Correct the background from the from factor
    FormFactor = S - (1-ModDepth)*B;
    FormFactor = FormFactor./(ModDepth*B);
    FormFactor = FormFactor/max(FormFactor);
    
    %Perform APT on background-corrected signal
    [FreqDimension,~] = size(K);
    FreqDistribution=zeros(1,FreqDimension);
    for k=1:FreqDimension % sum in eqn [21]
        FreqDistribution(k)=FreqDistribution(k)+sum(K(k,:).*FormFactor.*APT_t)/NormConstant(k);
    end
    APTdistribution = Crosstalk\FreqDistribution';
    
    %Get merit value for this background start value
    Merit(FitStartPos - StartPosMin + 1) = sum(abs(APTdistribution(1:3)));
end

[~,OptStartPos] = min(Merit);
FitStartPos = OptStartPos + StartPosMin - 1;
FitStartTime = t(FitStartPos);

%Store output result in the cache
Output = {FitStartTime,FitStartPos};
cachedData = addcache(cachedData,hashKey,Output);
