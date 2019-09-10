function [err,data,maxerr] = test(opt,oldata)


Ntime = 100;
Ndist = 200;

dt = 0.008;
t = linspace(0,dt*Ntime,Ntime);
[~,rmin,rmax] = time2dist(t);
r = linspace(rmin,rmax,Ndist);
InputParam = [3,0.5];
P = rd_onegaussian(r,[3,0.5]);

K = dipolarkernel(t,r);
DipEvoFcn = K*P;

InitialGuess = [2 0.1];
[FitParam,FitP] = fitparamodel(DipEvoFcn,@rd_onegaussian,r,K,InitialGuess);
err(1) = any(abs(FitP - P)>1e-5);
err(2) = any(abs(FitParam - InputParam)>1e-3);
err(3)  = length(FitP) < length(DipEvoFcn);
err = any(err);

maxerr = max(abs(FitP - P));
data = [];

if opt.Display
   figure(1),clf,hold on
   plot(t,DipEvoFcn,'b')
   plot(t,K*FitP,'r')
end

end