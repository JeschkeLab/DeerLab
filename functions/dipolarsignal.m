%
% DIPOLARSIGNAL Generate dipolar signal from distance distribution
%
%   [F,D] = DIPOLARSIGNAL(t,r,P)
%   Calculates the noiseless form factor (F) and dipolar evolution function (D)
%   on a time axis (t) from the distance distribution (P) in a distance axis (r).
%
%   [F,D] = DIPOLARSIGNAL(t,r)
%   If no distribution is provided, then the dipolar signal corresponding 
%   to a single distance (r) is computed.
%
%   [F,D] = DIPOLARSIGNAL(...,'Property',Value)
%   Additional options can be passed as property-value pairs. You can specify
%   several name and value pair arguments in any order.
%
%   Name-Value Pair Arguments:
%
%   'ModDepth' - Modulation depth of the form factor
%
%   'B' - Array containing a background function
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

function [FormFactor,S] = dipolarsignal(t,r,P,varargin)

if nargin<3
    P = [];
    if length(r)~=1
       error('Only one distance is accepted if P is not provided.') 
    end
end
if ischar(P)
    varargin = [{P} varargin];
    P = [];
end
%Parse optional input arguments
[ModDepth,B,NoiseLevel,Offset,Overtones] = parseoptional({'ModDepth','B','NoiseLevel','Offset','Overtones'},varargin);
%Validate inputs
if isempty(ModDepth)
    ModDepth = 1;
end
if isempty(B)
    B = ones(length(t),1);
end
if isempty(NoiseLevel)
    NoiseLevel = 0;
end
if isempty(Overtones)
    Overtones = 1;
end
if isempty(Offset)
    Offset = 1;
end
validateattributes(NoiseLevel,{'numeric'},{'scalar','nonnegative'},mfilename,'NoiseLevel')
validateattributes(ModDepth,{'numeric'},{'scalar','nonnegative','nonempty'},mfilename,'ModDepth')
validateattributes(t,{'numeric'},{'increasing','nonempty'},mfilename,'t')
validateattributes(r,{'numeric'},{'increasing','nonempty','nonnegative'},mfilename,'r')
validateattributes(B,{'numeric'},{'2d'},mfilename,'B')
validateattributes(P,{'numeric'},{'2d'},mfilename,'P')
validateattributes(Offset,{'numeric'},{'scalar','nonnegative'},mfilename,'Offset')
validateattributes(Overtones,{'numeric'},{'2d','nonnegative'},mfilename,'Overtones')

if ~isempty(P) && length(r)~=length(P)
    error('The distance axis and distribution lengths must be equal.')
end

if ModDepth>1 || ModDepth<0
    error('Modulation depth must be in the range of 0 to 1.')
end
if numel(unique(round(diff(r),12)))~=1 && length(r)~=1
    error('Distance axis must be a monotonically increasing vector.')
end


%Convert time step to microseconds if given in nanoseconds
usesNanoseconds = mean(diff(t))>=0.5;
if usesNanoseconds
    t = round(t)/1000; % ns->us
end

%Convert distance axis to nanoseconds if givne in Angstrom
if ~isnanometer(r)
    r = r/10;
end

%Get length of distribution
N = length(t);

%Normalize the distance distribution if not normalized

%Get the kernel
K = dipolarkernel(t,r,'OvertoneCoeffs',Overtones);

%Calculate dipolar evolution function
if ~isempty(P)
    P = P/sum(P)/mean(diff(r));
    S = K*P;
else
    S = K;
end

%Generate Gaussian noise
Noise = whitegaussnoise(N,NoiseLevel);

%Calculate form factor with backbround
FormFactor = (1-ModDepth) + ModDepth*S;
FormFactor = FormFactor.*B;

%Add noise and intensity offset
FormFactor = (FormFactor + Noise)*Offset;


end