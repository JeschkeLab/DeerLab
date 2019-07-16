function APTkernel = getAPTkernel(TimeAxis,varargin)

%Check if user requested some options via name-value input
[ExcitationBandwidth] = parseOptional({'ExcitationBandwidth'},varargin);
if ~isempty(ExcitationBandwidth)
    validateattributes(ExcitationBandwidth,{'numeric'},{'scalar','nonnegative'})
end
if iscolumn(TimeAxis)
   TimeAxis = TimeAxis'; 
end
validateattributes(TimeAxis,{'numeric'},{'nonempty','increasing','nonnegative'},'TimeAxis')

FreqElement = 1/(2*max(TimeAxis));
TimeDimension = length(TimeAxis);
FreqDimension = floor(TimeDimension/2)-2;
FreqAxis = linspace(1,FreqDimension,FreqDimension);
FreqAxis = FreqElement*(FreqAxis+1/4*ones(1,FreqDimension));

NormalizationFactor=zeros(1,FreqDimension); % initialize vector of normalization constant

%Numerical angular dipolar frequency
wdd=2*pi*FreqAxis';

%If given, account for limited excitation bandwidth
if ~isempty(ExcitationBandwidth)
    wdd = exp(-wdd.^2/ExcitationBandwidth^2).*wdd;
end

%Allocate products for speed
wddt = wdd.*TimeAxis;
kappa = sqrt(6*wddt/pi);

%Compute Fresnel integrals of 0th order
C = fresnelC(kappa);
S = fresnelS(kappa);

%Compute dipolar kernel
Base = sqrt(pi./(wddt*6)).*(cos(wddt).*C + sin(wddt).*S);
Base(:,1) = 1; 

%Normalize with respect to dipolar evolution time
for k=1:FreqDimension % normalize kernel traces to value at time origin
  Base(k,:) = Base(k,:)./Base(k,1);
  NormalizationFactor(k) = sum(Base(k,:).*Base(k,:).*TimeAxis); % compute normalization constant, eqn [19]
end

[FreqDimension,~] = size(Base); % size of kernel
Crosstalk = zeros(FreqDimension,FreqDimension); % initialize crosstalk matrix
for k=1:FreqDimension % compute crosstalk matrix, eqn [20]
  for l=1:FreqDimension
    mu = Base(k,:);
    Crosstalk(k,l) = sum(mu.*Base(l,:).*TimeAxis)/NormalizationFactor(k);
  end
end

%Construct the kernel object to be passed later to the APT.m function
APTkernel = aptkernel('Base',Base,...
                      'NormalizationFactor',NormalizationFactor,...
                      'FreqAxis',FreqAxis,...
                      'TimeAxis',TimeAxis,...
                      'Crosstalk',Crosstalk);
                  