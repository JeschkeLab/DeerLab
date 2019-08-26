function [err,data,maxerr] = test(opt,olddata)

clear regularize

%=======================================
% Check TV regularization
%=======================================
Dimension = 80;
TimeStep = 0.008;
TimeAxis = linspace(0,TimeStep*Dimension,Dimension);
DistanceAxis = time2dist(TimeAxis);
Distribution = onegaussian(DistanceAxis,[3,0.5]);

Kernel = dipolarkernel(TimeAxis,DistanceAxis);
DipEvoFcn = Kernel*Distribution;

%Set optimal regularization parameter (found numerically lambda=0.005)
RegParam = 0.0001;
RegMatrix = regoperator(Dimension,2);
Result = fitregmodel(DipEvoFcn,Kernel,DistanceAxis,RegMatrix,'huber',RegParam,'Solver','fmincon','HuberParam',1.35);

error = abs(Result - Distribution);
err(1) = any(error>2e-1);
maxerr = max(error);
err = any(err);
data = [];

if opt.Display
 	figure(8),clf
    hold on
    plot(DistanceAxis,Distribution,'k') 
    plot(DistanceAxis,Result,'r')
end

end