function [Base,NormalizationFactor,FreqAxis,TimeAxis,Crosstalk]=getAPTkernel(TimeDimension,TimeStep)
% Computes kernel for approximate Pake transformation
%
% numdat    number of data points in time domain
%
% base      kernel data
% norm      normalization constants, eqn [19] of Ref. (1)
% ny        dipolar frequency axis (MHz)
% t         time axis (µs)
% crosstalk cross-talk matrix, eqn [20] of Ref. 1
%
% (1) G. Jeschke, A. Koch, U. Jonas, A. Godt, J. Magn. Reson. 155, 72-82 (2002)
%
% (c) G. Jeschke, 2001,2019
%

if nargin<2 || isempty(TimeStep)
  TimeStep = 0.008;
end

TimeAxis = linspace(0,(TimeDimension-1)*TimeStep,TimeDimension);
FreqElement = 1/(2*max(TimeAxis));
FreqDimension = floor(TimeDimension/2)-2;
FreqAxis = linspace(1,FreqDimension,FreqDimension);
FreqAxis = FreqElement*(FreqAxis+1/4*ones(1,FreqDimension));

NormalizationFactor=zeros(1,FreqDimension); % initialize vector of normalization constant

wdd=2*pi*FreqAxis'; % angular frequency
%Allocate products for speed
wddt = wdd*TimeAxis;
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
