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