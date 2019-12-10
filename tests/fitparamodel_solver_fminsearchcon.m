function [err,data,maxerr] = test(opt,oldata)

rng(1)
Dimension = 200;
dt = 0.008;
t = linspace(0,dt*Dimension,Dimension);
r = time2dist(t);
paramin  = [3,0.5];
P = rd_onegaussian(r,paramin);

K = dipolarkernel(t,r);
S = K*P;

param0 = [2 0.1];
[fitparam,Pfit] = fitparamodel(S,@rd_onegaussian,r,K,param0,'Solver','fminsearchcon');
err(1) = any(abs(Pfit - P)>2e-5);
err(2) = any(abs(fitparam - paramin)>1e-3);
err = any(err);

maxerr = max(abs(Pfit - P));
data = [];

if opt.Display
   figure(1),clf,hold on
   plot(t,S,'b')
   plot(t,K*Pfit,'r')
end

end