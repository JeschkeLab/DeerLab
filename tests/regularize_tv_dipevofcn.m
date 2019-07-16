function [err,data] = test(opt,olddata)

%=======================================
% Check TV regularization
%=======================================
Dimension = 200;
TimeStep = 0.008;
TimeAxis = linspace(0,TimeStep*Dimension,Dimension);
DistanceAxis = time2distAxis(TimeAxis);
Distribution = gaussian(DistanceAxis,3,0.5);
Distribution = Distribution/sum(Distribution);

Kernel = getKernel(TimeAxis,DistanceAxis);
DipEvoFcn = Kernel*Distribution;

%Set optimal regularization parameter (found numerically lambda=0.005)
RegParam = 0.005;
TVResult1 = regularize(DipEvoFcn,Kernel,'tv',RegParam,'RegMatrixOrder',3,'Solver','fmincon');

err(1) = any(abs(TVResult1 - Distribution)>1e-2);
err = any(err);
data = [];

if opt.Display
 	figure(8),clf
    hold on
    plot(DistanceAxis,Distribution,'k') 
    plot(DistanceAxis,TVResult1,'r')
end

end