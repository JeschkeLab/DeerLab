%
% DIPOLARSIGNAL Generate dipolar signal from distance distribution
%
%   [F,D] = DIPOLARSIGNAL(t,r,P)
%   Calculates the noiseless form factor (F) and dipolar evolution function (D)
%   on a time axis (t) from the distance distribution (P) in a distance axis (r).
%
%   [F,D] = DIPOLARSIGNAL(...,'Property',Value)
%   Additional options can be passed as property-value pairs. You can specify
%   several name and value pair arguments in any order.
%
%   Name-Value Pair Arguments: 
%
%   'ModDepth' - Modulation depth of the form factor
%
%   'Background' - Array containing a background function
%
%   'NoiseLevel' - Level (standard deviation) of gaussian noise to add 
%
%   'Offset' - Vertical offset to add to the ouput signal
%
%   'Overtones' - Array of RIDME overtone coefficients 
%
%
% Copyright(C) 2019  Luis Fabregas, DeerAnalysis2
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License 3.0 as published by
% the Free Software Foundation.

function [FormFactor,DipEvoFcn] = dipolarsignal(TimeAxis,DistanceAxis,Distribution,varargin)

%Parse optional input arguments
[ModDepth,Background,NoiseLevel,Offset,Overtones] = parseoptional({'ModDepth','Background','NoiseLevel','Offset','Overtones'},varargin);
%Validate inputs
if isempty(ModDepth)
    ModDepth = 1;
end
if isempty(Background)
    Background = ones(length(TimeAxis),1);
end
if isempty(NoiseLevel)
    NoiseLevel = 0;
end
if isempty(Overtones)
    Overtones = 1;
end
if isempty(Offset)
    Offset = 0;
end
validateattributes(NoiseLevel,{'numeric'},{'scalar','nonnegative'},mfilename,'NoiseLevel')
validateattributes(ModDepth,{'numeric'},{'scalar','nonnegative','nonempty'},mfilename,'ModDepth')
validateattributes(TimeAxis,{'numeric'},{'increasing','nonempty'},mfilename,'TimeAxis')
validateattributes(DistanceAxis,{'numeric'},{'increasing','nonempty','nonnegative'},mfilename,'DistanceAxis')
validateattributes(Background,{'numeric'},{'2d'},mfilename,'Background')
validateattributes(Distribution,{'numeric'},{'2d','nonempty'},mfilename,'Background')
validateattributes(Offset,{'numeric'},{'scalar','nonnegative'},mfilename,'Offset')
validateattributes(Overtones,{'numeric'},{'2d','nonnegative'},mfilename,'Background')

if ModDepth>1 || ModDepth<0
    error('Modulation depth must be in the range of 0 to 1.')
end
if numel(unique(round(diff(DistanceAxis),12)))~=1
    error('Distance axis must be a monotonically increasing vector.')
end

%Get length of distribution
N = length(Distribution);

%Normalize the distance distribution if not normalized
Distribution = Distribution/sum(Distribution)/mean(diff(DistanceAxis));

%Get the kernel
Kernel = dipolarkernel(TimeAxis,DistanceAxis,[],'OvertoneCoeffs',Overtones);

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