function [err,data,maxerr] = test(opt,oldata)


Dimension = 200;
dt = 0.008;
t = linspace(0,dt*Dimension,Dimension);
r = time2dist(t);
InputParam = [3 0.5];
P = rd_onegaussian(r,[3,0.5]);

K = dipolarkernel(t,r);
DipEvoFcn = K*P;

InitialGuess = [2 0.1];
[FitParam,FitP] = fitparamodel(DipEvoFcn,@rd_onegaussian,r,K,'Solver','fmincon');
err(1) = any(abs(FitP - P)>1e-5);
err(2) = any(abs(FitParam - InputParam)>1e-3);
err = any(err);

maxerr = max(abs(FitP - P));
data = [];

if opt.Display
   figure(1),clf,hold on
   plot(t,DipEvoFcn,'b')
   plot(t,K*FitP,'r')
end

end