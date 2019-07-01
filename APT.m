function [APTDistr] = APT(Signal,Opts)

if ~isa(Signal,'pdsdata')
  error('First argument must a valid pdsdata class')
end
if nargin>1
  if ~isa(Opts,'daopts')
    error('First argument must a valid daopts class')
  end
end
if isempty(Signal.ExpData)
  error('ExpData property is empty.')
end
if isempty(Signal.TimeAxis)
  error('TimeAxis property is empty.')
end


[Kernel,NormConstant,APT_FrequencyAxis,APT_TimeAxis,Crosstalk] = getAPTkernel(Signal.Length,Signal.TimeStep/1000);

DDS = Opts.DistDomainSmoothing;

[FreqDimension,~] = size(Kernel); % size of kernel
FreqDistribution = zeros(1,FreqDimension); % initialize distribution
for k=1:FreqDimension % sum in eqn [21]
    FreqDistribution(k) = FreqDistribution(k)+sum(Kernel(k,:).*Signal.DipEvoFcn.*APT_TimeAxis)/NormConstant(k);
end
APTdistribution = Crosstalk\FreqDistribution'; % crosstalk correction, eqn [22]

Freq2Dist = zeros(length(APT_FrequencyAxis),1); % initialize distance axis (mapping of dipolar frequencies to distances)
for k = 1:length(Freq2Dist)
    Freq2Dist(k) = (52.04/APT_FrequencyAxis(k))^(1/3);
end

Freq2DistExt = APT_FrequencyAxis/2; % initialize distance axis (mapping of dipolar frequencies to distances)
for k = 1:length(Freq2DistExt)
    Freq2DistExt(k) = (52.04/Freq2DistExt(k))^(1/3);
end

APTdistribution = interp1(Freq2Dist,APTdistribution,Freq2DistExt,'pchip',0);

APTdistribution = APTdistribution';

FilteredAPTdistribution = zeros(length(APTdistribution),1);
for k = 1:length(APTdistribution) % smoothing, convolute with...
    DDSfilter = (Freq2DistExt-Freq2DistExt(k)*ones(1,length(Freq2DistExt)))/DDS; % ... Gaussian line as filter
    DDSfilter = exp(-DDSfilter.*DDSfilter);
    DDSfilter = DDSfilter';
    FilteredAPTdistribution(k) = sum(DDSfilter.*APTdistribution)/sum(DDSfilter); % normalization!
end

for k=1:length(APTdistribution)
    FilteredAPTdistribution(k) = FilteredAPTdistribution(k)/(Freq2DistExt(k))^4; % keep constant integral after mapping to distances
end

FilteredAPTdistribution = FilteredAPTdistribution/max(FilteredAPTdistribution); % renormalize distribution

FitDipEvolFcn = Kernel'*FilteredAPTdistribution;
FitDipEvolFcn = FitDipEvolFcn/FitDipEvolFcn(1);

APTDistr = distdistr;
APTDistr.Distribution = FilteredAPTdistribution;
APTDistr.DistanceAxis = Freq2DistExt;
APTDistr.ExpSignal = Signal.DipEvoFcn;
APTDistr.FitSignal = FitDipEvolFcn;
APTDistr.TimeAxis = Signal.TimeAxis;

% 
% APTDistribution = get_std_distr(Freq2DistExt,FilteredAPTdistribution,APT_DistanceAxis);
% APTDistribution = 0.01*APTDistribution/sum(APTDistribution);
% 
% tpcf=handles.Pake_t;
% 
% distr0=get_std_distr(Freq2DistExt/stretch,FilteredAPTdistribution,APT_DistanceAxis);
% distr0=0.01*distr0/sum(distr0);
% deer=pcf2deer(distr0,kernel,APT_DistanceAxis,tpcf);
% % bsl=sum(deer(925:1024))/100;
% % deer=deer-bsl*ones(size(deer));
% deer=deer/max(deer);
% sim=interp1(tpcf*TimeStep/0.008,deer,tdip);