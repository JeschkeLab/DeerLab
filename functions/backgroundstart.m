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
%   'SearchStart' - Time [us] at which the background start
%                   search starts (default=0.1*max(t)).
%
%   'SearchEnd' - Time [us] at which the backgrund start search
%                 stops (default=0.6*max(t)).
%
%   'EndCutOff' - Time [us] after which the signal is no longer used for fitting
%

% This file is a part of DeerLab. License is MIT (see LICENSE.md). 
% Copyright(c) 2019: Luis Fabregas, Stefan Stoll, Gunnar Jeschke and other contributors.

function [FitStartTime,FitStartPos] = backgroundstart(V,t,BckgModel,varargin)


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

if ~iscolumn(t)
    t = t';
end
if ~iscolumn(V)
    V = V.';
end
if ~isreal(V)
    error('Input signal cannot be complex.')
end
validateattributes(V,{'numeric'},{'2d','nonempty'},mfilename,'FitData')
validateattributes(t,{'numeric'},{'2d','nonempty','increasing'},mfilename,'t')

%--------------------------------------------------------------------------
% Parse & Validate Optional Input
%--------------------------------------------------------------------------
%Check if user requested some options via name-value input
[SearchStart,SearchEnd,EndCutoff] = parseoptional({'RelSearchStart','RelSearchEnd','EndCutoffPos'},varargin);

if isempty(SearchStart)
    SearchStart = 0.1*max(t);
end

if isempty(SearchEnd)
    SearchEnd = 0.6*max(t);
end
if  isempty(EndCutoff)
    EndCutoff = max(t);
else
    validateattributes(EndCutoff,{'numeric'},{'scalar','nonempty'},mfilename,'EndCutoffPos')
end

if SearchEnd > max(t)
    error('SearchEnd option exceeds the largest value in the time-axis.')
end

if SearchStart < min(t)
    error('SearchStart option is smaller than the smallest value in the time-axis.')
end

if SearchEnd<SearchStart
   error('SearchStart option cannot be larger than SearchEnd option.') 
end

%--------------------------------------------------------------------------
% Adaptive background correction start search
%--------------------------------------------------------------------------
[~,EndCutoffPos] = min(abs(t - EndCutoff));
t = t(1:EndCutoffPos);
V = V(1:EndCutoffPos);

%Get APT kernel
APTkernel = aptkernel(t);

%Get APT kernel data
K = APTkernel.Base;
NormConstant = APTkernel.NormalizationFactor;
APT_t = APTkernel.t(:).';
Crosstalk = APTkernel.Crosstalk;

%Search for background fit start
[~,StartPosMin] = min(abs(t - SearchStart));
[~,StartPosMax] = min(abs(t - SearchEnd));
if StartPosMax > EndCutoffPos
   StartPosMax = EndCutoffPos; 
end

%Preallocate Merit variable
Merit = zeros(1,StartPosMax-StartPosMin);

for FitStartPos = StartPosMin:StartPosMax
    
    %Define data to be fitted according to current background start
    FitStart = t(FitStartPos);
    
    %Fit the background with current start
    [B,lambda] = fitbackground(V,t,BckgModel,FitStart);

    %Correct the background from the from factor
    F = V - (1-lambda)*B;
    F = F./(lambda*B);
    F = F/max(F);
    
    %Perform APT on background-corrected signal
    [FreqDimension,~] = size(K);
    FreqP=zeros(1,FreqDimension);
    for k=1:FreqDimension % sum in eqn [21]
        FreqP(k)=FreqP(k)+sum(K(k,:).*F.'.*APT_t)/NormConstant(k);
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
