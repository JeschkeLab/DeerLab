%
% TIME2FREQ Conversion from time-axis to frequency-axis
%
%   nu = TIME2FREQ(t)
%   Computes the N-point frequency axis (nu) from the N-point input 
%   time axis (t) according to the Nyquist criterion.
%
%   nu = TIME2FREQ(t,M)
%   Computes the (M)-point frequency axis (nu) from the N-point input 
%   time axis (t) according to the Nyquist criterion.
%
% Copyright(C) 2019  Luis Fabregas, DeerAnalysis2
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License 3.0 as published by
% the Free Software Foundation.

function FrequencyAxis = time2freq(TimeAxis,FreqPoints)

if nargin<2 
    FreqPoints = length(TimeAxis);
end

TimeStep = mean(diff(TimeAxis));
FrequencyAxis = linspace(-1/(2*TimeStep),1/(2*TimeStep),FreqPoints);

end