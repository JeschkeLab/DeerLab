function [err,data,maxerr] = test(opt,olddata)

%=======================================
% Check TV regularization
%=======================================
Dimension = 200;
TimeStep = 0.008;
TimeAxis = linspace(0,TimeStep*Dimension,Dimension);
DistanceAxis = time2dist(TimeAxis);
Distribution = rd_onegaussian(DistanceAxis,[3,0.5]);

Kernel = dipolarkernel(TimeAxis,DistanceAxis);
DipEvoFcn = Kernel*Distribution;

%Set optimal regularization parameter (found numerically lambda=0.005)
RegParam = 0.0005;
RegMatrix = regoperator(Dimension,3);
TVResult1 = fitregmodel(DipEvoFcn,Kernel,DistanceAxis,RegMatrix,'tv',RegParam,'Solver','fmincon');

error = abs(TVResult1 - Distribution);
err(1) = any(error>2e-1);
maxerr = max(error);
err = any(err);
data = [];

if opt.Display
 	figure(8),clf
    hold on
    plot(DistanceAxis,Distribution,'k') 
    plot(DistanceAxis,TVResult1,'r')
end

end