function [err,data,maxerr] = test(opt,oldata)


warning('off','all')

Dimension = 200;
dt = 0.008;
t = linspace(0,dt*Dimension,Dimension);
r = time2dist(t);
InputParam = [3 0.5];
P = rd_onegaussian(r,[3,0.5]);

K = dipolarkernel(t,r);
S = K*P;

mymodel = @(t,param)K*rd_onegaussian(r,param);

InitialGuess = [2 0.1];
[Sfit,FitParam] = fitparamodel(S,mymodel,t,InitialGuess);
Pfit = rd_onegaussian(r,FitParam);
err(1) = any(abs(Pfit - P)>1e-5);
err(2) = any(abs(FitParam - InputParam)>1e-3);
err = any(err);

maxerr = max(abs(Pfit - P));
data = [];

warning('on','all')

if opt.Display
   figure(1),clf,hold on
   plot(t,S,'b')
   plot(t,Sfit,'r')
end

end