function [err,data,maxerr] = test(opt,oldata)


Dimension = 200;
TimeStep = 0.008;
TimeAxis = linspace(0,TimeStep*Dimension,Dimension);
DistanceAxis = time2dist(TimeAxis);
InputParam1 = [3 0.5];
InputParam2 = [4 0.5];
mixedModel = mixmodels({@onegaussian,@onegaussian});
mixedmodelParameters = [0.3 InputParam1 InputParam2];
MixedDistribution = mixedModel(DistanceAxis,mixedmodelParameters);

Kernel = dipolarkernel(TimeAxis,DistanceAxis);
Signal = Kernel*MixedDistribution;

Fit = fitparamodel(Signal,Kernel,DistanceAxis,mixedModel,[]);

err = any(abs(MixedDistribution - Fit)>1e-5);
maxerr = max(abs(MixedDistribution - Fit));
data = [];

if opt.Display
   figure(1),clf,hold on
   plot(TimeAxis,MixedDistribution,'b')
   plot(TimeAxis,Fit,'r')
end

end