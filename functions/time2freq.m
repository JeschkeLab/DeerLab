function FrequencyAxis = time2freq(TimeAxis,FreqPoints)

if nargin<2 
    FreqPoints = length(Signal);
end

TimeStep = mean(diff(TimeAxis));
FrequencyAxis = linspace(-1/(2*TimeStep),1/(2*TimeStep),FreqPoints);

end