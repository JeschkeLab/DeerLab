function data = prepareSignal(data,options)
% Perfmorm basic pre-processing of PDS primary data: phase correction,
% zero-time correction, background fitting, and construction of the
% form factor and dipolar evolution function.

if ~isa(data,'pdsdata')
  error('First argument must a valid pdsdata class')
end
if nargin>1
  if ~isa(options,'optda')
    error('First argument must a valid optda class')
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
[data.ClusterFcn,data.Phase] = correctPhase(data.ClusterFcn);
[data.ClusterFcn,data.TimeAxis,data.ZeroTime] = correctZeroTime(data.ClusterFcn,data.TimeAxis);

%Fit background
[FitStartTime,FitStartPos] = getBackgroundStart(data.ClusterFcn ,data.TimeAxis);
data.FitStart = FitStartTime;
Data2fit = data.ClusterFcn(FitStartPos:end);
FitTimeAxis = data.TimeAxis(FitStartPos:end);

data.Background = fitBackground(Data2fit,data.TimeAxis,FitTimeAxis,'exponential',[]);

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