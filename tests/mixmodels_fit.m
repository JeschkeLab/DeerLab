function [err,data,maxerr] = test(opt,oldata)


Dimension = 200;
dt = 0.008;
t = linspace(0,dt*Dimension,Dimension);
r = time2dist(t);
InputParam1 = [3 0.5];
InputParam2 = [4 0.5];
mixedModel = mixmodels({@rd_onegaussian,@rd_onegaussian});
mixedmodelParameters = [0.3 InputParam1 InputParam2];
MixedP = mixedModel(r,mixedmodelParameters);

K = dipolarkernel(t,r);
S = K*MixedP;

[~,Fit] = fitparamodel(S,mixedModel,r,K);

err = any(abs(MixedP - Fit)>1e-5);
maxerr = max(abs(MixedP - Fit));
data = [];

if opt.Display
   figure(1),clf,hold on
   plot(t,MixedP,'b')
   plot(t,Fit,'r')
end

end