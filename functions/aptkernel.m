% 
% APTKERNEL Computes the dipolar interaction kernel and elements required 
%           for the approximate Pake transformation (APT).  
%
%       K = APTKERNEL(t) 
%       Computes a structure (K) containing the (N/2-2)xN point kernel,
%       the (N/2-2) point array of normalization factors, N/2-2) point
%       frequency axis and the (N/2-2)x(N/2-2) crosstalk matrix corresponding
%       to the N-point time axis (t).
%
%       K = APTKERNEL(T,'ExcitationBandwidth',w)
%       The excitation bandwidth (w) of the experiment can be passed as an
%       option to account for it in the kernel.
%

% This file is a part of DeerAnalysis. License is MIT (see LICENSE.md). 
% Copyright(c) 2019: Luis Fabregas, Stefan Stoll, Gunnar Jeschke and other contributors.


function APTkernel = aptkernel(t,varargin)


%Check if user requested some options via name-value input
[ExcitationBandwidth] = parseoptional({'ExcitationBandwidth'},varargin);
if ~isempty(ExcitationBandwidth)
    validateattributes(ExcitationBandwidth,{'numeric'},{'scalar','nonnegative'})
end
if iscolumn(t)
   t = t'; 
end
validateattributes(t,{'numeric'},{'nonempty','increasing'},'t')

%Convert time step to microseconds if given in nanoseconds
usesNanoseconds = mean(diff(t))>=0.5;
if usesNanoseconds
    t = t/1000; % ns->us
end
%Use absolute time scale, required for negative times
t = abs(t);

%--------------------------------------------------------------------------
%Memoization
%--------------------------------------------------------------------------

persistent cachedData
if isempty(cachedData)
    cachedData =  java.util.LinkedHashMap;
end
hashKey = datahash({t,varargin});
if cachedData.containsKey(hashKey)
    Output = cachedData.get(hashKey);
    [APTkernel] = java2mat(Output);
    APTkernel.NormalizationFactor = APTkernel.NormalizationFactor.';
    APTkernel.FreqAxis = APTkernel.FreqAxis.';
    APTkernel.t = APTkernel.t.';
    return
end

%--------------------------------------------------------------------------
%------------------------------------------------------------------------


%Turn off warnings to avoid ill-conditioned warnings 
warning('off','all')

FreqElement = 1/(2*max(t));
TimeDimension = length(t);
FreqDimension = floor(TimeDimension/2)-2;
FreqAxis = linspace(1,FreqDimension,FreqDimension);
FreqAxis = FreqElement*(FreqAxis+1/4*ones(1,FreqDimension));

NormalizationFactor = zeros(1,FreqDimension); % initialize vector of normalization constant

%Numerical angular dipolar frequency
wdd = 2*pi*FreqAxis';

%Allocate products for speed
wddt = wdd.*t;
kappa = sqrt(6*wddt/pi);

%Compute Fresnel integrals of 0th order
C = fresnelC(kappa);
S = fresnelS(kappa);

%Compute dipolar kernel
Base = sqrt(pi./(wddt*6)).*(cos(wddt).*C + sin(wddt).*S);
Base(isnan(Base)) = 1; 

%If given, account for limited excitation bandwidth
if ~isempty(ExcitationBandwidth)
    Base = exp(-wdd'.^2/ExcitationBandwidth^2)'.*Base;
end

Base = Base*mean(abs(diff(FreqAxis)));

%Normalize with respect to dipolar evolution time
for k=1:FreqDimension % normalize kernel traces to value at time origin
  NormalizationFactor(k) = sum(Base(k,:).*Base(k,:).*t); % compute normalization constant, eqn [19]
end

[FreqDimension,~] = size(Base); % size of kernel
Crosstalk = zeros(FreqDimension,FreqDimension); % initialize crosstalk matrix
for k=1:FreqDimension % compute crosstalk matrix, eqn [20]
  for l=1:FreqDimension
    mu = Base(k,:);
    Crosstalk(k,l) = sum(mu.*Base(l,:).*t)/NormalizationFactor(k);
  end
end

%Construct the kernel object to be passed later to the APT.m function
APTkernel = struct('Base',Base,...
                      'NormalizationFactor',NormalizationFactor,...
                      'FreqAxis',FreqAxis,...
                      't',t,...
                      'Crosstalk',Crosstalk);

%Store output result in the cache
cachedData = addcache(cachedData,hashKey,APTkernel);

%Turn warnings back on
warning('on','all')

end          