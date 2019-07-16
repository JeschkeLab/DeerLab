function [DistanceAxis,rmin,rmax] = time2distAxis(TimeAxis)

validateattributes(TimeAxis,{'numeric'},{'nonempty','increasing','nonnegative'})

TimeStep = mean(diff(TimeAxis));
Dimension = length(TimeAxis);
ny0 = 52.04;
rmin = (4*TimeStep*ny0/0.85)^(1/3);
rmax = 6*(Dimension*TimeStep/2)^(1/3);

DistanceAxis = linspace(rmin,rmax,Dimension);

end