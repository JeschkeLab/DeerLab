function [err,data,maxerr] = test(opt,olddata)

%=======================================
% Check Tikhonov regularization
%=======================================

Dimension = 100;
TimeStep = 0.008;
TimeAxis = linspace(0,TimeStep*Dimension,Dimension);
DistanceAxis = time2dist(TimeAxis);
Distribution = gaussian(DistanceAxis,3,0.5);
Distribution = Distribution/sum(Distribution);

Kernel = dipolarkernel(TimeAxis,DistanceAxis);
DipEvoFcn = Kernel*Distribution;

%Set optimal regularization parameter (found numerically lambda=0.13)
RegMatrix = regoperator(Dimension,2);
RegParamRange = regparamrange(Kernel,RegMatrix);
% RegParam = selregparam(RegParamRange,DipEvoFcn,Kernel,RegMatrix,'gcv','regtype','tikhonov','NonNegConstrained',false);
RegParam = 100;
Result = regularize(DipEvoFcn,Kernel,RegMatrix,'tikhonov',RegParam,'NonNegConstrained',false);

err(1) = any(abs(Result - Distribution)>1e-2);
maxerr = max(abs(Result - Distribution));
err = any(err);
data = [];

if opt.Display
 	figure(8),clf
    hold on
    plot(DistanceAxis,Distribution,'k') 
    plot(DistanceAxis,Result,'r')
end

end