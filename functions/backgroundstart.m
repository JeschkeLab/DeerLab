function [FitStartTime,FitStartPos] = backgroundstart(Signal,TimeAxis,EndCutoffPos,BckgModel,ModelParam,varargin)


%--------------------------------------------------------------------------
% Parse & Validate Required Input
%--------------------------------------------------------------------------
if nargin<2
    error('Not enough input arguments.')
end

if nargin<3 || isempty(EndCutoffPos)
    EndCutoffPos = length(TimeAxis);
else
    validateattributes(EndCutoffPos,{'numeric'},{'scalar','nonempty'},mfilename,'EndCutoffPos')
end

if nargin<3 || nargin<4 || isempty(BckgModel)
    BckgModel = 'exponential';
    ModelParam = [];
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
[RelSearchStart,RelSearchEnd] = parseoptional({'RelSearchStart','RelSearchEnd'},varargin);

if isempty(RelSearchStart)
    RelSearchStart = 0.1;
else
    
end

if isempty(RelSearchEnd)
    RelSearchEnd = 0.6;
else
    
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
