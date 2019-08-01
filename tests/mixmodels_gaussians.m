function [err,data,maxerr] = test(opt,oldata)


Dimension = 200;
TimeStep = 0.008;
TimeAxis = linspace(0,TimeStep*Dimension,Dimension);
DistanceAxis = time2dist(TimeAxis);
InputParam1 = [3 0.5];
Distribution1 = onegaussian(DistanceAxis,InputParam1);
Distribution1 = Distribution1/sum(Distribution1);
InputParam2 = [4 0.5];
Distribution2 = onegaussian(DistanceAxis,InputParam2);
Distribution2 = Distribution2/sum(Distribution2);

Distribution = 0.7*Distribution2 + 0.3*Distribution1;

mixedModel = mixmodels({@onegaussian,@onegaussian});

mixedmodelParameters = [0.3 InputParam1 InputParam2];

MixedDistribution = mixedModel(DistanceAxis,mixedmodelParameters);

err = any(abs(MixedDistribution - Distribution)>1e-8);
maxerr = max(abs(MixedDistribution - Distribution));
data = [];

if opt.Display
   figure(1),clf,hold on
   plot(TimeAxis,Distribution,'b')
   plot(TimeAxis,MixedDistribution,'r')
end

end