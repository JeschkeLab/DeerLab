function [err,data] = test(opt,olddata)

%=======================================
% Test supressGhostDistances.m
%=======================================


Dimension = 100;
TimeAxis = linspace(0.006,3,Dimension);
DistanceAxis = time2dist(TimeAxis);

Distribution = gaussian(DistanceAxis,3,0.5);
DistrAB = Distribution/sum(Distribution)';
Distribution = gaussian(DistanceAxis,4,0.5);
DistrAC = Distribution/sum(Distribution)';
Distribution = gaussian(DistanceAxis,5,0.5);
DistrBC = Distribution/sum(Distribution)';
Distribution = DistrAB + DistrBC + DistrAC;
Distribution = Distribution/sum(Distribution)';

FreqAxis = 52.04./(DistanceAxis.^3)';

wddt = 2*pi*FreqAxis.*TimeAxis;
kappa = sqrt(6*wddt/pi);
%Compute Fresnel integrals of 0th order
C = fresnelC(kappa);
S = fresnelS(kappa);

omegaAB = sum(sqrt(pi./(wddt*6)).*(cos(wddt).*C + sin(wddt).*S).*DistrAB);
omegaBC = sum(sqrt(pi./(wddt*6)).*(cos(wddt).*C + sin(wddt).*S).*DistrBC);
omegaAC = sum(sqrt(pi./(wddt*6)).*(cos(wddt).*C + sin(wddt).*S).*DistrAC);

lambda = 0.6;

P = 1/3*(omegaAB + omegaBC + omegaAC);
T = 1/6*(omegaAB.*omegaAC + omegaAB.*omegaBC + omegaBC.*omegaAC);


FormFactor3 = 1-2*lambda + lambda^2 +lambda*(1-lambda)*P + lambda^2*T;
FormFactor2 = 1-lambda + lambda*P;

correctedFormFactor3 = FormFactor3.^(1/2);

Dip3 = (correctedFormFactor3-(1-lambda))/lambda;
Dip2 = (FormFactor2-(1-lambda))/lambda;

Dip2 = Dip2/Dip2(1);
Dip3 = Dip3/Dip3(1);

signal = supressGhostDistances(FormFactor3,lambda,3);


Kernel = getKernel(TimeAxis,DistanceAxis);
RegParam = 4;
Dist2 = regularize(Dip2,Kernel,'tikhonov',RegParam,'Solver','fnnls');
Dist3 = regularize(Dip3,Kernel,'tikhonov',RegParam,'Solver','fnnls');
DistrTest = regularize(signal,Kernel,'tikhonov',RegParam,'Solver','fnnls');

err(1) = any(abs(DistrTest - Dist3)>1e-6);
err(2) = any(abs(DistrTest - Dist2)>1e-2);

err = any(err);
data = [];


end

