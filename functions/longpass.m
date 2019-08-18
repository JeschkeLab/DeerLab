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
% Copyright(C) 2019  Luis Fabregas, DeerAnalysis2
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License 3.0 as published by
% the Free Software Foundation.

function FilteredSignal = longpass(TimeAxis,Signal,PassDist,Steepness)

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
validateattributes(TimeAxis,{'numeric'},{'nonempty','increasing'},mfilename,'TimeAxis')
validateattributes(Signal,{'numeric'},{'nonempty'},mfilename,'Signal')
if iscolumn(Signal)
    Signal = Signal';
end


%Get the absolute time axes and timestep
TimeStep = mean(abs(diff(TimeAxis)));

%Compute width of transition band in distance domain as in MATLAB's lowpass function
TransitionBand = (1 - Steepness)*(PassDist);
StopDist = PassDist - 1/2*TransitionBand;

%Get all frequencies in MHz
PassBandFreq = 52.04/(PassDist^3);
StopBandFreq = 52.04/(StopDist^3);
SamplingFreq = 1/(TimeStep);

%Convert all frequencies to Hz
PassBandFreq = PassBandFreq*1e6;
StopBandFreq = StopBandFreq*1e6;
SamplingFreq = SamplingFreq*1e6;

%Filter the signal using a windowed lowpass FIR filter
FilteredSignal = winlowpass(Signal,StopBandFreq,PassBandFreq,SamplingFreq);

end
