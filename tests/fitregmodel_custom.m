function [err,data,maxerr] = test(opt,olddata)

%=======================================
% Check Tikhonov regularization
%=======================================

Dimension = 200;
TimeStep = 0.008;
TimeAxis = linspace(0,TimeStep*Dimension,Dimension);
DistanceAxis = time2dist(TimeAxis);
Distribution = gaussian(DistanceAxis,3,0.5);
Distribution = Distribution/sum(Distribution)/mean(diff(DistanceAxis));

Kernel = dipolarkernel(TimeAxis,DistanceAxis);
DipEvoFcn = Kernel*Distribution;


%Set optimal regularization parameter (found numerically lambda=0.13)
RegParam = 0.02;
RegMatrix = regoperator(Dimension,3);
RegFunctional = @(Dist)(1/2*norm(Kernel*Dist - DipEvoFcn)^2 + RegParam^2*max(RegMatrix*Dist)^2);
Result = fitregmodel(DipEvoFcn,Kernel,DistanceAxis,RegMatrix,RegFunctional,RegParam,'Solver','fmincon');

err = any(abs(Result - Distribution)>1e-1);
maxerr = max(abs(Result - Distribution));

data = [];

if opt.Display
 	figure(8),clf
    hold on
    plot(DistanceAxis,Distribution,'k') 
    plot(DistanceAxis,Result,'r')
end

end