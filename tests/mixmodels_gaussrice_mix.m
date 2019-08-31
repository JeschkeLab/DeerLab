function [err,data,maxerr] = test(opt,oldata)


Dimension = 200;
dt = 0.008;
t = linspace(0,dt*Dimension,Dimension);
r = time2dist(t);
InputParam1 = [2.5 0.3];
Distribution1 = rd_onegaussian(r,InputParam1);
Distribution1 = Distribution1/sum(Distribution1)/mean(diff(r));
InputParam2 = [3.5 0.3];
Distribution2 = rd_onegaussian(r,InputParam2);
Distribution2 = Distribution2/sum(Distribution2)/mean(diff(r));

InputParam3 = [4.5 0.3];
Distribution3 = rd_onerice(r,InputParam3);
Distribution3 = Distribution3/sum(Distribution3)/mean(diff(r));

Distribution = 0.4*Distribution2 + 0.3*Distribution1 + 0.3*Distribution3;

mixedModel = mixmodels({@rd_onegaussian,@rd_onegaussian,@rd_onerice});

mixedmodelParameters = [0.3 0.4 InputParam1 InputParam2 InputParam3];

MixedDistribution = mixedModel(r,mixedmodelParameters);

err = any(abs(MixedDistribution - Distribution)>1e-8);
maxerr = max(abs(MixedDistribution - Distribution));
data = [];

if opt.Display
   figure(1),clf,hold on
   plot(t,Distribution,'b')
   plot(t,MixedDistribution,'r')
end

end