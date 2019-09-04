function [err,data,maxerr] = test(opt,oldata)


Dimension = 200;
dt = 0.008;
t = linspace(0,dt*Dimension,Dimension);
r = time2dist(t);
InputParam = [3 0.4 4.5 0.3 0.6];
P = rd_tworice(r,InputParam);
P = P/sum(P)/mean(diff(r));

K = dipolarkernel(t,r);
DipEvoFcn = K*P;

InitialGuess = [2 0.1 1 0.6 0.2];
[FitP] = fitparamodel(DipEvoFcn,@rd_tworice,r,K,InitialGuess,'solver','fmincon');
err = any(abs(FitP - P)>1e-5);

maxerr = max(abs(FitP - P));
data = [];

if opt.Display
   figure(1),clf
   subplot(121)
   hold on
   plot(t,DipEvoFcn,'b')
   plot(t,K*FitP,'r')
   subplot(122)
   hold on
   plot(r,P,'b')
   plot(r,FitP,'r')
end

end