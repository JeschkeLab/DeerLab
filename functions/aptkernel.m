% 
% APTKERNEL Computes the dipolar interaction kernel and elements required 
%           for the approximate Pake transformation (APT).  
%
%       K = APTKERNEL(T) 
%       Computes a structure K containing the (N/2-2)xN point kernel,
%       the (N/2-2) point array of normalization factors, N/2-2) point
%       frequency axis and the (N/2-2)x(N/2-2) crosstalk matrix corresponding
%       to the N-point time axis T.
%
%       K = APTKERNEL(T,'ExitationBandwidth',W)
%       The excitation bandwidth of the experiment can be passed as an
%       option to account for it in the kernel.
%
% Copyright(C) 2019  Luis Fabregas, DeerAnalysis2
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License 3.0 as published by
% the Free Software Foundation.

function APTkernel = aptkernel(TimeAxis,varargin)

%Check if user requested some options via name-value input
[ExcitationBandwidth] = parseoptional({'ExcitationBandwidth'},varargin);
if ~isempty(ExcitationBandwidth)
    validateattributes(ExcitationBandwidth,{'numeric'},{'scalar','nonnegative'})
end
if iscolumn(TimeAxis)
   TimeAxis = TimeAxis'; 
end
validateattributes(TimeAxis,{'numeric'},{'nonempty','increasing'},'TimeAxis')
TimeAxis = abs(TimeAxis);

%--------------------------------------------------------------------------
%Memoization
%--------------------------------------------------------------------------

persistent cachedData
if isempty(cachedData)
    cachedData =  java.util.LinkedHashMap;
end
hashKey = datahash({TimeAxis,varargin});
if cachedData.containsKey(hashKey)
    Output = cachedData.get(hashKey);
    [APTkernel] = java2mat(Output);
    APTkernel.NormalizationFactor = APTkernel.NormalizationFactor.';
    APTkernel.FreqAxis = APTkernel.FreqAxis.';
    APTkernel.TimeAxis = APTkernel.TimeAxis.';
    return
end

%--------------------------------------------------------------------------
%------------------------------------------------------------------------

FreqElement = 1/(2*max(TimeAxis));
TimeDimension = length(TimeAxis);
FreqDimension = floor(TimeDimension/2)-2;
FreqAxis = linspace(1,FreqDimension,FreqDimension);
FreqAxis = FreqElement*(FreqAxis+1/4*ones(1,FreqDimension));

NormalizationFactor = zeros(1,FreqDimension); % initialize vector of normalization constant

%Numerical angular dipolar frequency
wdd = 2*pi*FreqAxis';

%Allocate products for speed
wddt = wdd.*TimeAxis;
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
APTkernel = struct('Base',Base,...
                      'NormalizationFactor',NormalizationFactor,...
                      'FreqAxis',FreqAxis,...
                      'TimeAxis',TimeAxis,...
                      'Crosstalk',Crosstalk);

%Store output result in the cache
cachedData = addcache(cachedData,hashKey,APTkernel);


end          