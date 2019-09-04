function [err,data,maxerr] = test(opt,oldata)


Dimension = 200;
dt = 0.008;
t = linspace(0,dt*Dimension,Dimension);
r = time2dist(t);
InputParam1 = [3 0.5];
P1 = rd_onegaussian(r,InputParam1);
P1 = P1/sum(P1)/mean(diff(r));
InputParam2 = [4 0.5];
P2 = rd_onegaussian(r,InputParam2);
P2 = P2/sum(P2)/mean(diff(r));

P = 0.7*P2 + 0.3*P1;

mixedModel = mixmodels({@rd_onegaussian,@rd_onegaussian});

mixedmodelParameters = [0.3 InputParam1 InputParam2];

MixedP = mixedModel(r,mixedmodelParameters);

err = any(abs(MixedP - P)>1e-8);
maxerr = max(abs(MixedP - P));
data = [];

if opt.Display
   figure(1),clf,hold on
   plot(t,P,'b')
   plot(t,MixedP,'r')
end

end