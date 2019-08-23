%
% TIME2DIST Conversion from time-axis to distance-axis
%
%   [r,rmin,rmax] = TIME2DIST(t)
%   Computes the N-point distance axis (r) according to the input time axis (t).
%   The minimal and maximal distances (rmin,rmax) are determined by the empirical
%   approximations derived by Gunnar Jeschke as implemented in the older
%   DeerAnalysis versions.
%
%   [r,rmin,rmax] = TIME2DIST(t,M)
%   The length of the output axis can be specified by the additional input
%   parameter (M).
%
% Copyright(C) 2019  Luis Fabregas, DeerAnalysis2
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License 3.0 as published by
% the Free Software Foundation.

function [DistanceAxis,rmin,rmax] = time2dist(TimeAxis,M)

validateattributes(TimeAxis,{'numeric'},{'nonempty','increasing'})

if nargin<2 || isempty(M)
    M = length(TimeAxis);
else
   validateattributes(M,{'numeric'},{'scalar','nonnegative'}) 
end

TimeAxis = abs(TimeAxis);
TimeStep = mean(abs(diff(TimeAxis)));
tmax = max(TimeAxis);
ny0 = 52.04;
rmin = (4*TimeStep*ny0/0.85)^(1/3);
rmax = 6*(tmax/2)^(1/3);

DistanceAxis = linspace(rmin,rmax,M);

end