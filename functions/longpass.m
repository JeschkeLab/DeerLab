function FilteredSignal = longpass(TimeAxis,Signal,PassBandDist)
% Low-pass filtering to suppress proton modulations, suppresses distances
% below threshold r_filter_min (default: 1.5 nm)

if nargin < 3 || isempty(PassBandDist)
    PassBandDist = 1.5;
else
validateattributes(PassBandDist,{'numeric'},{'nonempty','scalar','nonnegative'},mfilename,'PassBandDist')
end

validateattributes(TimeAxis,{'numeric'},{'nonempty','increasing','nonnegative'},mfilename,'TimeAxis')
validateattributes(Signal,{'numeric'},{'nonempty'},mfilename,'Signal')

if iscolumn(Signal)
    Signal = Signal';
end

%Get frequencies in MHz
PassBandFreq = 52.04/(PassBandDist^3);
TimeStep = mean(diff(TimeAxis));
SamplingFreq = 1/TimeStep;

%Get frequencies in Hz
PassBandFreq = PassBandFreq*1e6;
SamplingFreq = SamplingFreq*1e6;

%Apply low-pass filter
FilteredSignal = lowpass(Signal,PassBandFreq,SamplingFreq);


end
