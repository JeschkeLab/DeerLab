function [pass,maxerr] = test(opt)


Dimension = 200;
dt = 0.008;
t = linspace(0,dt*Dimension,Dimension);
r = time2dist(t);
InputParam1 = [2.5 0.3];
P1 = rd_onegaussian(r,InputParam1);
P1 = P1/sum(P1)/mean(diff(r));
InputParam2 = [3.5 0.3];
P2 = rd_onegaussian(r,InputParam2);
P2 = P2/sum(P2)/mean(diff(r));

InputParam3 = [4.5 0.3];
P3 = rd_onerice(r,InputParam3);
P3 = P3/sum(P3)/mean(diff(r));

P = 0.4*P2 + 0.3*P1 + 0.3*P3;

mixedModel = mixmodels({@rd_onegaussian,@rd_onegaussian,@rd_onerice});

mixedmodelParameters = [0.3 0.4 InputParam1 InputParam2 InputParam3];

MixedP = mixedModel(r,mixedmodelParameters);

pass = all(abs(MixedP - P)>1e-8);
maxerr = max(abs(MixedP - P));
 

if opt.Display
   figure(1),clf,hold on
   plot(t,P,'b')
   plot(t,MixedP,'r')
end

end