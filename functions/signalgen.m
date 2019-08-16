%
% SIGNALGEN Generate signal from distance distribution
%
%   [F,D] = SIGNALGEN(T,R,P,B,LAMBDA)
%   Calculates the noiseless form factor F and dipolar evolution function D
%   from the distance distribution P in a distance axis R. The background B
%   must be defined on a time axis T and the modulation depth can be passed
%   as LAMBDA.
%
%   [F,D] = SIGNALGEN(T,R,P,B,LAMBDA,NOISELEVEL)
%   Adds gaussian noise to the calculated signals with a noise standard 
%   deviation given by the NOISELEVEL input argument.
%
%   [F,D] = SIGNALGEN(T,R,P,B,LAMBDA,NOISELEVEL,OFFSET)
%   Adds and additional vertical offset to the signal given by OFFSET
%   
%
% Copyright(C) 2019  Luis Fabregas, DeerAnalysis2
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License 3.0 as published by
% the Free Software Foundation.

function [FormFactor,DipEvoFcn] = signalgen(TimeAxis,DistanceAxis,Distribution,Background,ModDepth,NoiseLevel,Offset)


if nargin<6 || isempty(NoiseLevel)
    NoiseLevel = 0;
end

if nargin<7 || isempty(Offset)
    Offset = 0;
end

%Validate input
validateattributes(NoiseLevel,{'numeric'},{'scalar','nonnegative'},mfilename,'NoiseLevel')
validateattributes(ModDepth,{'numeric'},{'scalar','nonnegative','nonempty'},mfilename,'ModDepth')
validateattributes(TimeAxis,{'numeric'},{'increasing','nonempty'},mfilename,'TimeAxis')
validateattributes(DistanceAxis,{'numeric'},{'increasing','nonempty','nonnegative'},mfilename,'DistanceAxis')
validateattributes(Background,{'numeric'},{'2d'},mfilename,'Background')
validateattributes(Distribution,{'numeric'},{'2d','nonempty'},mfilename,'Background')
validateattributes(Offset,{'numeric'},{'scalar','nonnegative'},mfilename,'Offset')

%Get length of distribution
N = length(Distribution);

if isempty(Background)
    Background = ones(N,1);
end

if ModDepth>1 || ModDepth<0 
   error('Modulation depth must be in the range of 0 to 1.')
end
if numel(unique(round(diff(DistanceAxis),12)))~=1
    error('Distance axis must be a monotonically increasing vector.')
end

%Normalize the distance distribution if not normalized
Distribution = Distribution/sum(Distribution)/mean(diff(DistanceAxis));

%Get the kernel
Kernel = dipolarkernel(TimeAxis,DistanceAxis);

%Calculate dipolar evolution function
DipEvoFcn = Kernel*Distribution;

%Generate Gaussian noise
Noise = whitenoise(N,NoiseLevel);

%Calculate form factor with backbround
FormFactor = (1-ModDepth) + ModDepth*DipEvoFcn;
FormFactor = FormFactor.*Background;

%Add noise and intensity offset
FormFactor = FormFactor + Noise + Offset;


end