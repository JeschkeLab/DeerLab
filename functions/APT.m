function [APTDistribution] = APT(DipEvoFcn,APTkernel,DistDomainSmoothing)

if ~isa(APTkernel,'aptkernel')
    error('The input APTkernel must be a a valid aptkernel class object.')
end
if iscolumn(DipEvoFcn)
   DipEvoFcn = DipEvoFcn'; 
end

%Get APT kernel data
[Kernel,NormConstant,APT_FrequencyAxis,APT_TimeAxis,Crosstalk] = dismountAPTkernel(APTkernel);

%Compute frequency distribution
[FreqDimension,~] = size(Kernel);
FreqDistribution = zeros(1,FreqDimension);
for k=1:FreqDimension
    FreqDistribution(k) = FreqDistribution(k)+sum(Kernel(k,:).*DipEvoFcn.*APT_TimeAxis)/NormConstant(k);
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
    DDSfilter = (MappedDistances - MappedDistances(k)*ones(1,length(MappedDistances)))/DistDomainSmoothing;
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
UniformDistanceAxis = linspace(min(MappedDistances),max(MappedDistances),length(DipEvoFcn));
APTDistribution = uniformGrain(MappedDistances,FilteredAPTdistribution,UniformDistanceAxis);


APTDistribution = apt('DistanceAxis',UniformDistanceAxis,...
                      'Distribution',APTDistribution,...
                      'DistDomainSmoothing',DistDomainSmoothing);
end

