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

% This file is a part of DeerAnalysis. License is MIT (see LICENSE.md). 
% Copyright(c) 2019: Luis Fabregas, Stefan Stoll, Gunnar Jeschke and other contributors.

function [FitStartTime,FitStartPos] = backgroundstart(S,t,BckgModel,varargin)


%Turn off warnings to avoid ill-conditioned warnings at each iteration
warning('off','all')


%--------------------------------------------------------------------------
% Parse & Validate Required Input
%--------------------------------------------------------------------------
if nargin<3
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
if ~isreal(S)
    error('Input signal cannot be complex.')
end
validateattributes(S,{'numeric'},{'2d','nonempty'},mfilename,'FitData')
validateattributes(t,{'numeric'},{'2d','nonempty','increasing'},mfilename,'t')

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

if RelSearchEnd<RelSearchStart
   error('RelSearchStart option cannot be larger than RelSearchEnd option.') 
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
    FreqP=zeros(1,FreqDimension);
    for k=1:FreqDimension % sum in eqn [21]
        FreqP(k)=FreqP(k)+sum(K(k,:).*FormFactor.*APT_t)/NormConstant(k);
    end
    APTdistribution = Crosstalk\FreqP';
    
    %Get merit value for this background start value
    Merit(FitStartPos - StartPosMin + 1) = sum(abs(APTdistribution(1:3)));
end

[~,OptStartPos] = min(Merit);
FitStartPos = OptStartPos + StartPosMin - 1;
FitStartTime = t(FitStartPos);

%Turn warnings back on
warning('on','all')
