function [FilteredSignal,FilteredSpectrum] = longpassFilter(TimeAxis,Signal,PassBandDist)
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

temp = Signal;

%Adjust length of singal and axis to even number
if mod(length(Signal),2) == 1
    Signal = Signal(1:length(Signal) - 1);
    TimeAxis = TimeAxis(1:length(TimeAxis) - 1);
end

%Pass-band frequency in MHz
PassBandFreq = 52.04/(PassBandDist^3);

%Construct frequency axis
TimeStep = mean(diff(TimeAxis));
FreqAxis = linspace(-1/(2*TimeStep),1/(2*TimeStep),length(TimeAxis));
	
%Get start/end indices of the filter
[~,FilterStart] = min(abs(FreqAxis + PassBandFreq));
[~,FilterEnd] = min(abs(FreqAxis - PassBandFreq));
	
%Get spectrum
Spectrum=fftshift(fft(Signal));
	
FilteredSpectrum = Spectrum;
FitData = FilteredSpectrum(1:FilterStart - 1);
FitAxis = linspace(1,length(FitData),length(FitData));
[Fit,~] = polyfit(FitAxis,FitData,3);
Correction = polyval(Fit,FitAxis);

FilteredSpectrum(1:FilterStart - 1) = Correction;

FitData = FilteredSpectrum(FilterEnd:length(FilteredSpectrum));
FitAxis = linspace(1,length(FitData),length(FitData));
[Fit,~] = polyfit(FitAxis,FitData,3);
Correction = polyval(Fit,FitAxis);

FilteredSpectrum(FilterEnd:length(FilteredSpectrum)) = Correction;

FilteredSignal = ifft(fftshift(FilteredSpectrum));
FilteredSignal = real(FilteredSignal);	

if length(temp)>length(FilteredSignal)
    FilteredSignal=[FilteredSignal temp(length(temp))];
end


end
