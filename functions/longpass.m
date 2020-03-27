%
% LONGPASS Filter short distances by means of a lowpass FIR filter
%
%   X = LONGPASS(t,S)
%   Applies a lowpass filter to the signal (S) with a N-point time axis (t) to 
%   suppress any distances below 1.5nm. 
%
%   X = LONGPASS(t,S,rpass)
%   Applies a lowpass filter to the signal (S) with a N-point time axis (t) to 
%   suppress any distances below a given pass distance (rpass) in nm. 
%
%   X = LONGPASS(t,S,rpass,st)
%   Applies a lowpass filter to the signal (S) with a N-point time axis (t) to 
%   suppress any distances below a given pass distance (rpass) in nm. The
%   steepness of the transition band in distance domain can be passed as (st).
%

% This file is a part of DeerLab. License is MIT (see LICENSE.md). 
% Copyright(c) 2019: Luis Fabregas, Stefan Stoll, Gunnar Jeschke and other contributors.

function FilteredS = longpass(t,S,PassDist,Steepness)

if nargin < 3 || isempty(PassDist)
    PassDist = 1.5;
end
if nargin<4 || isempty(Steepness)
    Steepness  = 0.80;
end
if Steepness>1
    error('Steepness must be a scalar value between 0 and 1.')
end
validateattributes(Steepness,{'numeric'},{'nonempty','scalar','nonnegative'},mfilename,'Steepness')
validateattributes(PassDist,{'numeric'},{'nonempty','scalar','nonnegative'},mfilename,'PassBandDist')
validateattributes(t,{'numeric'},{'nonempty','increasing'},mfilename,'t')
validateattributes(S,{'numeric'},{'nonempty'},mfilename,'S')
S = S(:);
if ~isreal(S)
    error('Input signal cannot be complex.')
end

%Get the absolute time axes and timestep
dt = mean(abs(diff(t)));

%Compute width of transition band in distance domain as in MATLAB's lowpass function
TransitionBand = (1 - Steepness)*(PassDist);
StopDist = PassDist - 1/2*TransitionBand;

%Get all frequencies in MHz
PassBandFreq = 52.04/(PassDist^3);
StopBandFreq = 52.04/(StopDist^3);
SamplingFreq = 1/dt;

%Convert all frequencies to Hz
PassBandFreq = PassBandFreq*1e6;
StopBandFreq = StopBandFreq*1e6;
SamplingFreq = SamplingFreq*1e6;


%Filter the signal using a windowed lowpass FIR filter
FilteredS = winlowpass(S,StopBandFreq,PassBandFreq,SamplingFreq,'ForwardBackward',true);


end
