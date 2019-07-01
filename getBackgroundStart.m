function [FitStartTime,FitStartPos]=getBackgroundStart(Signal,TimeAxis,EndCutoffPos,BckgModel,ModelParam)

if nargin<5 || isempty(BckgModel)
  BckgModel = 'exponential';
  ModelParam = [];
end

if nargin<4 || isempty(EndCutoffPos)
  EndCutoffPos = length(TimeAxis);
end

%Get zero-time position
[~,ZeroTimePosition] = min(abs(TimeAxis));
TimeAxis = TimeAxis(ZeroTimePosition:EndCutoffPos);
TimeStep = TimeAxis(2) - TimeAxis(1);
Signal = real(Signal((ZeroTimePosition:EndCutoffPos)));
Length = length(TimeAxis);

[Kernel,NormConstant,~,APT_TimeAxis,Crosstalk] = getAPTkernel(Length,);

% imagflag=get(handles.imaginary,'Value');
% dcmplx=imagflag*handles.cmplx; % Complex data?
%
% % reference deconvolution for variable-time DEER, normalization for
% % constant-time DEER
% if handles.ctvt,
%     z=handles.A_vexp;
%     ref=Signal(1,:);
%     sc=max(real(ref));
%     sig=Signal(2,:);
%     Signal=real(sig)./real(ref)+1i*imag(ref)/sc;
% else
% 	Signal=Signal/max(real(Signal));
% end;

% Adaptive background correction
%Search for background fit start between 10% and 60% of the signal
StartPosMin = round(0.1*Length);
if StartPosMin<1
  StartPosMin=1;
end
StartPosMax = round(0.6*Length);
if StartPosMax<5
  StartPosMax=5;
end

Merit = zeros(1,StartPosMax-StartPosMin);

for FitStartPos=StartPosMin:StartPosMax
  
  FitTimeAxis=TimeAxis(FitStartPos:length(TimeAxis)); % time window of baseline region
  FitData=Signal(FitStartPos:length(Signal)); % experimental data in this window
  
  % Background fit
  Background=fitBackground(FitData,TimeAxis,FitTimeAxis,BckgModel,ModelParam);
  
  FormFactor=Signal-Background; % subtract background
  FormFactor=FormFactor./Background; % divide by background, eqn [13]
  
  FormFactor = FormFactor/max(FormFactor); % normalize
  [FreqDimension,~] = size(Kernel); % size of kernel
  FreqDistribution=zeros(1,FreqDimension); % initialize distribution
  for k=1:FreqDimension % sum in eqn [21]
    FreqDistribution(k)=FreqDistribution(k)+sum(Kernel(k,:).*FormFactor.*APT_TimeAxis)/NormConstant(k);
  end
  APTdistribution = Crosstalk\FreqDistribution'; % crosstalk correction, eqn [22]
  Merit(FitStartPos - StartPosMin + 1) = sum(abs(APTdistribution(1:3)));
end

[~,OptStartPos] = min(Merit);
FitStartPos = OptStartPos + StartPosMin - 1;
FitStartTime = TimeAxis(FitStartPos);
