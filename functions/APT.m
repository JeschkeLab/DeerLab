function [APTDistr] = APT(Signal,Opts)

if ~isa(Signal,'DAsignal')
  error('First argument must a valid DAsignal class object.')
end
if nargin>1
  if ~isa(Opts,'DAoptions')
    error('Second argument must a valid DAoptions class object.')
  end
end
if isempty(Signal.DipEvoFcn)
  error('Signal has not been prepared.')
end
if isempty(Signal.TimeAxis)
  error('TimeAxis property is empty.')
end

%Get APT kernel data
[Kernel,NormConstant,APT_FrequencyAxis,APT_TimeAxis,Crosstalk] = getAPTkernel(Signal.Length,Signal.TimeStep/1000);

%Compute frequency distribution
[FreqDimension,~] = size(Kernel);
FreqDistribution = zeros(1,FreqDimension);
for k=1:FreqDimension
    FreqDistribution(k) = FreqDistribution(k)+sum(Kernel(k,:).*Signal.DipEvoFcn.*APT_TimeAxis)/NormConstant(k);
end

%Perform crosstalk correction
APTdistribution = Crosstalk\FreqDistribution'; % crosstalk correction, eqn [22]

%Map fequencies to distances
Freq2Dist = zeros(length(APT_FrequencyAxis),1); % initialize distance axis (mapping of dipolar frequencies to distances)
for k = 1:length(Freq2Dist)
    Freq2Dist(k) = (52.04/APT_FrequencyAxis(k))^(1/3);
end
MappedDistances = APT_FrequencyAxis/2; % initialize distance axis (mapping of dipolar frequencies to distances)
for k = 1:length(MappedDistances)
    MappedDistances(k) = (52.04/MappedDistances(k))^(1/3);
end
APTdistribution = interp1(Freq2Dist,APTdistribution,MappedDistances,'pchip',0);
APTdistribution = APTdistribution';

%Perform distance-domain smoothing filtering
FilteredAPTdistribution = zeros(length(APTdistribution),1);
for k = 1:length(APTdistribution)
    DDSfilter = (MappedDistances - MappedDistances(k)*ones(1,length(MappedDistances)))/Opts.DistDomainSmoothing;
    DDSfilter = exp(-DDSfilter.*DDSfilter);
    DDSfilter = DDSfilter';
    FilteredAPTdistribution(k) = sum(DDSfilter.*APTdistribution)/sum(DDSfilter);
end

%Normalize integral
for k=1:length(APTdistribution)
    FilteredAPTdistribution(k) = FilteredAPTdistribution(k)/(MappedDistances(k))^4;
end
%Renormalize
FilteredAPTdistribution = FilteredAPTdistribution/max(FilteredAPTdistribution); 

%Interpolate to uniform distance axis
UniformDistanceAxis = linspace(min(MappedDistances),max(MappedDistances),length(Signal.TimeAxis));                 
APTDistribution = uniformGrain(MappedDistances,FilteredAPTdistribution,UniformDistanceAxis);

%Construct distance distribution class output
APTDistr = DAdistribution('DistanceAxis',UniformDistanceAxis,...
                          'Distribution',APTDistribution,...
                          'Signal',Signal.DipEvoFcn,...
                          'TimeAxis',Signal.TimeAxis,...
                          'SignalID',Signal.ID);
end

