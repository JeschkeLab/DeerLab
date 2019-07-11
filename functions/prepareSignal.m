function data = prepareSignal(data,opts)
% Perfmorm basic pre-processing of PDS primary data: phase correction,
% zero-time correction, background fitting, and construction of the
% form factor and dipolar evolution function.

if ~isa(data,'DAsignal')
  error('First argument must a valid DAsignal class object')
end
if nargin>1
  if ~isa(opts,'DAoptions')
    error('First argument must a valid DAoptions class object')
  end
end
if isempty(data.ExpData)
  error('ExpData property is empty.')
end
if isempty(data.TimeAxis)
  error('TimeAxis property is empty.')
end

%Normalize cluster signal
data.ClusterFcn = data.ExpData/data.ExpData(1);
[data.ClusterFcn,data.Phase] = correctPhase(data.ClusterFcn,opts.Phase);
[data.ClusterFcn,data.TimeAxis,data.ZeroTime] = correctZeroTime(data.ClusterFcn,data.TimeAxis,opts.ZeroTime);

%Fit background
if isempty(opts.Background)
if strcmp(opts.FitStartSearch,'auto')
    [FitStartTime,FitStartPos] = getBackgroundStart(data.ClusterFcn ,data.TimeAxis,opts.EndCutoffPos,opts.BackgroundModel,opts.ModelParam);
else
    FitStartPos = opts.FitStartPos;
    FitStartTime = data.TimeAxis(FitStartPos);
end
data.FitStart = FitStartTime;
Data2fit = data.ClusterFcn(FitStartPos:end);
FitTimeAxis = data.TimeAxis(FitStartPos:end);

data.Background = fitBackground(Data2fit,data.TimeAxis,FitTimeAxis,opts.BackgroundModel,opts.ModelParam);
else
   data.Background =  opts.Background;
end
%Correct for background by division
data.FormFactor = data.ClusterFcn./data.Background;
data.FormFactor = data.FormFactor/data.FormFactor(1);

%Calculate modulation depth
data.ModDepth = 1 - data.Background(1);

%Get dipolar evoution function
DipolarEvolution = data.FormFactor - (1 - data.ModDepth);
DipolarEvolution = DipolarEvolution./DipolarEvolution(1);
data.DipEvoFcn = DipolarEvolution;
end