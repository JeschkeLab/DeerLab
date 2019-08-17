%
% TIME2DIST Conversion from time-axis to distance-axis
%
%   [r,rmin,rmax] = TIME2DIST(t)
%   Computes the distance axis (r) according to the input time axis (t).
%   The minimal and maximal distances (rmin,rmax) are determined by the empirical
%   approximations derived by Gunnar Jeschke as implemented in the older
%   DeerAnalysis versions.
%
% Copyright(C) 2019  Luis Fabregas, DeerAnalysis2
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License 3.0 as published by
% the Free Software Foundation.

function [DistanceAxis,rmin,rmax] = time2dist(TimeAxis)

validateattributes(TimeAxis,{'numeric'},{'nonempty','increasing'})
TimeAxis = abs(TimeAxis);
TimeStep = mean(abs(diff(TimeAxis)));
tmax = max(TimeAxis);
ny0 = 52.04;
rmin = (4*TimeStep*ny0/0.85)^(1/3);
rmax = 6*(tmax/2)^(1/3);

DistanceAxis = linspace(rmin,rmax,length(TimeAxis));

end