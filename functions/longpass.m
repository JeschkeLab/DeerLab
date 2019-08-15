function FilteredSignal = longpass(TimeAxis,Signal,PassDist)
% Low-pass filtering to suppress proton modulations, suppresses distances
% below threshold r_filter_min (default: 1.5 nm)

if nargin < 3 || isempty(PassDist)
    PassDist = 1.5;
else
validateattributes(PassDist,{'numeric'},{'nonempty','scalar','nonnegative'},mfilename,'PassBandDist')
end

validateattributes(TimeAxis,{'numeric'},{'nonempty','increasing'},mfilename,'TimeAxis')
validateattributes(Signal,{'numeric'},{'nonempty'},mfilename,'Signal')

if iscolumn(Signal)
    Signal = Signal';
end
Steepness  = 0.80;
TimeAxis = abs(TimeAxis);
TimeStep = mean(abs(diff(TimeAxis)));

%Compute width of transition band as in MATLAB's lowpass function
TransitionBand = (1 - Steepness)*(PassDist);
StopDist = PassDist - 1/2*TransitionBand;

%Get frequencies in MHz
PassBandFreq = 52.04/(PassDist^3);
StopBandFreq = 52.04/(StopDist^3);
TimeStep = mean(abs(diff(TimeAxis)));
SamplingFreq = 1/(TimeStep);

%Get frequencies in Hz
PassBandFreq = PassBandFreq*1e6;
StopBandFreq = StopBandFreq*1e6;
SamplingFreq = SamplingFreq*1e6;

 FilteredSignal = winlowpass(Signal,StopBandFreq,PassBandFreq,SamplingFreq);



end
