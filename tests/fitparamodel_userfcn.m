function [err,data,maxerr] = test(opt,oldata)

warning('off','all')

Ntime = 200;
Ndist = 200;

dt = 0.008;
t = linspace(0,dt*Ntime,Ntime);
[~,rmin,rmax] = time2dist(t);
r = linspace(rmin,rmax,Ndist);
InputParam = [3,0.2];
P = rd_onegaussian(r,[3,0.2]);

K = dipolarkernel(t,r);
S = K*P;

InitialGuess = [3.5 0.3];
fcnhandle = @(r,param)exp(-((r-param(1))/(param(2))).^2);

[FitParam,FitP] = fitparamodel(S,fcnhandle,r,K,InitialGuess);
err(1) = any(abs(FitP - P)>1e-5);
err(2)  = length(FitP) < length(S);
err = any(err);

maxerr = max(abs(FitP - P));
data = [];

warning('on','all')

if opt.Display
   figure(1),clf,hold on
   plot(t,P,'b')
   plot(t,FitP,'r')
end

end