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

% This file is a part of DeerLab. License is MIT (see LICENSE.md). 
% Copyright(c) 2019: Luis Fabregas, Stefan Stoll, Gunnar Jeschke and other contributors.


function nu = time2freq(t,FreqPoints)

if nargin<2 
    FreqPoints = length(t);
end

dt = mean(diff(t));
nu = linspace(-1/(2*dt),1/(2*dt),FreqPoints);
nu = nu(:);

end